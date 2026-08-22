//! Reading and writing multi-element configurations.
//!
//! Two formats, with two different jobs:
//!
//! - **Multi-block `.dat`** — the XFOIL / MSES coordinate format, as used by
//!   `flexfoil-ui/public/airfoils/30p-30n.dat`. It carries element
//!   *coordinates* only: no placement, no paneling parameters, no reference
//!   quantities. [`Configuration::from_dat_blocks`] builds a configuration from
//!   coordinate blocks already split out of such a file;
//!   [`Configuration::to_dat`] writes one back out, separating elements with
//!   the `999.0 999.0` sentinel.
//! - **JSON** — the configuration format proper,
//!   [`Configuration::to_json`] / [`Configuration::from_json`]. It carries the
//!   geometry *and* the placements, the per-element paneling parameters and the
//!   reference quantities, so a configuration is a saveable artefact rather
//!   than something reassembled by hand from separate coordinate files.
//!
//! # Splitting a `.dat` file into blocks happens elsewhere
//! Deciding where one element's coordinates end and the next one's begin is a
//! lexical question about a file — blank lines, comment lines, `999.0`
//! sentinels, and stray fragments that belong to their neighbour. That parser
//! lives in `rustfoil-cli` (`load_airfoil_blocks`), with a matching one in
//! `packages/flexfoil-python` and one in `flexfoil-ui`. This module consumes
//! the blocks it produces rather than re-deriving them, so there is one set of
//! separator rules per language and not two.
//!
//! Those three parsers agreeing by inspection rather than by construction is a
//! known seam: a single splitter in this crate, with the CLI and the Python and
//! TypeScript bindings calling into it, would remove it. That is a follow-up,
//! not part of this module.
//!
//! # Why the JSON is written by hand
//! This crate depends only on `nalgebra`, deliberately — it compiles to WASM
//! and the bundle size is a product constraint. None of the types a
//! configuration is made of derive `serde::Serialize` / `Deserialize`, and
//! [`Point`] is an alias for `nalgebra::Point2<f64>`, whose serde support is a
//! `nalgebra` feature. Deriving would therefore mean a new dependency plus
//! edits to `body.rs`, `placement.rs`, `spline.rs` and `configuration.rs`. The
//! reader and writer here are self-contained instead, and the field layout is
//! the one a derive would produce, so adding derives later is a
//! format-compatible change.
//!
//! # Numeric fidelity
//! Both writers format every coordinate with `{:?}`, which is the shortest
//! decimal that reads back as the same `f64`. A configuration therefore
//! survives a write/read cycle with every coordinate bit-identical, in either
//! format — no `{:.6}`-style truncation of imported geometry.

use crate::body::Body;
use crate::configuration::{Configuration, Element};
use crate::error::GeometryError;
use crate::panel::Panel;
use crate::placement::Placement;
use crate::point::{point, vec2, Point};
use crate::spline::PanelingParams;

/// XFOIL / MSES element separator sentinel.
///
/// A coordinate line whose **both** values reach this is an element boundary,
/// never a point — the same rule the `.dat` parsers apply on the way in.
pub const DAT_SEPARATOR_SENTINEL: f64 = 999.0;

/// The element separator line [`Configuration::to_dat`] writes.
pub const DAT_SEPARATOR_LINE: &str = " 999.0 999.0";

/// The `format` field [`Configuration::to_json`] writes, and the only value
/// [`Configuration::from_json`] accepts.
pub const JSON_FORMAT_TAG: &str = "flexfoil.configuration";

/// The `version` field [`Configuration::to_json`] writes.
///
/// [`Configuration::from_json`] reads this version and earlier, and reports a
/// newer one rather than guessing at fields it does not know.
pub const JSON_FORMAT_VERSION: u32 = 1;

/// Role name given to the largest-chord element of an imported multi-element
/// file — see [`Configuration::from_dat_blocks`].
const ROLE_MAIN: &str = "main";

/// Role name given to elements forward of the main element.
const ROLE_SLAT: &str = "slat";

/// Role name given to elements aft of the main element.
const ROLE_FLAP: &str = "flap";

/// The contour a [`Body`] was built from.
///
/// `Body` keeps its geometry as panels, so the node list is recovered as every
/// panel's start point plus the last panel's end point. That is exactly the
/// slice [`Body::from_points`] was given, including the duplicated final point
/// of a closed (sharp trailing edge) contour, so
/// `Body::from_points(&body.name, &body_contour(&body))` reproduces the body.
///
/// Empty for a body with no panels.
///
/// # Example
/// ```
/// use rustfoil_core::body::Body;
/// use rustfoil_core::config_io::body_contour;
/// use rustfoil_core::point::point;
///
/// let points = vec![
///     point(1.0, 0.0),
///     point(0.5, -0.05),
///     point(0.0, 0.0),
///     point(0.5, 0.05),
///     point(1.0, 0.0),
/// ];
/// let body = Body::from_points("diamond", &points).unwrap();
/// assert_eq!(body_contour(&body), points);
/// ```
pub fn body_contour(body: &Body) -> Vec<Point> {
    let panels: &[Panel] = body.panels();
    match panels.last() {
        None => Vec::new(),
        Some(last) => {
            let mut contour = Vec::with_capacity(panels.len() + 1);
            contour.extend(panels.iter().map(|panel| panel.p1));
            contour.push(last.p2);
            contour
        }
    }
}

impl Configuration {
    /// Build a configuration from the per-element coordinate blocks of a
    /// multi-block `.dat` file.
    ///
    /// `name` is the file's name line; `blocks` are its element coordinate
    /// blocks in file order, as produced by the `.dat` parsers described in the
    /// [module documentation](crate::config_io). One [`Element`] is built per
    /// block, in the same order, each with an identity [`Placement`] and
    /// default [`PanelingParams`] — the format carries neither.
    ///
    /// # Element names
    /// A `.dat` file states no roles, so they are derived the same geometric
    /// way decision D4 picks the main element:
    ///
    /// - One block: identical to [`Configuration::single`] — the element takes
    ///   its name from the body, which takes it from `name`.
    /// - Several blocks: the largest-chord element is `"main"`, elements whose
    ///   chord midpoint lies forward of the main element's are `"slat"`, the
    ///   rest are `"flap"`. Where more than one element shares a role they are
    ///   numbered in file order (`"slat 1"`, `"slat 2"`). Body names are
    ///   `"<name> element <n>"`, `n` counting from 1, since every block came
    ///   from the same source file.
    ///
    /// These are defaults, not deductions. A caller that knows the roles should
    /// overwrite [`Element::name`].
    ///
    /// # Errors
    /// - [`GeometryError::ParseError`] if `blocks` is empty.
    /// - Whatever [`Body::from_points`] reports for a block that is not a
    ///   usable contour — too few points, or a repeated point giving a
    ///   zero-length panel.
    ///
    /// # Example
    /// ```
    /// use rustfoil_core::configuration::Configuration;
    /// use rustfoil_core::point::point;
    ///
    /// // Two blocks, as a `.dat` parser would hand them over.
    /// let contour = |x0: f64, chord: f64| vec![
    ///     point(x0 + chord, 0.0),
    ///     point(x0 + 0.5 * chord, -0.05 * chord),
    ///     point(x0, 0.0),
    ///     point(x0 + 0.5 * chord, 0.05 * chord),
    ///     point(x0 + chord, 0.0),
    /// ];
    /// let blocks = vec![contour(0.0, 1.0), contour(0.95, 0.3)];
    ///
    /// let config = Configuration::from_dat_blocks("two element", &blocks).unwrap();
    ///
    /// assert_eq!(config.len(), 2);
    /// assert_eq!(config.element(0).unwrap().name, "main");
    /// assert_eq!(config.element(1).unwrap().name, "flap");
    /// assert!(config.element(1).unwrap().placement.is_identity());
    /// ```
    pub fn from_dat_blocks(name: &str, blocks: &[Vec<Point>]) -> Result<Self, GeometryError> {
        if blocks.is_empty() {
            return Err(GeometryError::ParseError {
                line: 0,
                message: "no coordinate blocks: a configuration needs at least one element",
            });
        }

        let single = blocks.len() == 1;
        let mut elements = Vec::with_capacity(blocks.len());
        for (index, block) in blocks.iter().enumerate() {
            let body_name = if single {
                name.to_string()
            } else {
                format!("{name} element {}", index + 1)
            };
            elements.push(Element::from_body(Body::from_points(&body_name, block)?));
        }

        let mut config = Self::new(elements);
        if !single {
            assign_default_roles(&mut config);
        }
        Ok(config)
    }

    /// Write the configuration as a multi-block `.dat` file.
    ///
    /// `name` becomes the file's name line. Elements are written in
    /// configuration order, separated by [`DAT_SEPARATOR_LINE`], each
    /// coordinate at full `f64` precision (see the [module
    /// documentation](crate::config_io)).
    ///
    /// # Placements are baked in
    /// The `.dat` format has nowhere to put a placement, so each element's
    /// coordinates are written **with its placement applied** — the file
    /// describes where the elements actually are, which is the only thing it
    /// can describe faithfully. Reading it back therefore gives the same
    /// geometry carried by identity placements rather than the original
    /// placement decomposition. Use [`to_json`](Configuration::to_json) to keep
    /// the placements. For a configuration whose placements are already the
    /// identity — anything imported from a `.dat` file — this distinction does
    /// not arise, and the coordinates are passed through bit for bit.
    ///
    /// # The name line
    /// A `.dat` reader takes the first line as a name only when it does not
    /// parse as a coordinate pair. `name` is written verbatim when that holds,
    /// so a single-element file keeps its usual `NACA 0012` header; otherwise
    /// (an empty name, or one that reads as two numbers) it is written as a
    /// `#` comment so that it cannot be mistaken for geometry. Any embedded
    /// newline is replaced by a space, since the name occupies one line.
    ///
    /// # Errors
    /// [`GeometryError::InvalidParameter`] for a coordinate that cannot be
    /// written faithfully: one that is not finite, or one whose x *and* y both
    /// reach [`DAT_SEPARATOR_SENTINEL`] and would read back as an element
    /// separator.
    pub fn to_dat(&self, name: &str) -> Result<String, GeometryError> {
        let mut out = String::new();
        out.push_str(&dat_name_line(name));
        out.push('\n');

        for (index, element) in self.elements.iter().enumerate() {
            if index > 0 {
                out.push_str(DAT_SEPARATOR_LINE);
                out.push('\n');
            }
            for node in body_contour(&element.body) {
                out.push_str(&dat_coordinate_line(element.placement.apply(node))?);
                out.push('\n');
            }
        }

        Ok(out)
    }

    /// Write the configuration as JSON: geometry, placements, per-element
    /// paneling parameters and reference quantities.
    ///
    /// The document is a stable, reviewable shape — one coordinate pair per
    /// line — tagged with [`JSON_FORMAT_TAG`] and [`JSON_FORMAT_VERSION`]:
    ///
    /// ```text
    /// {
    ///   "format": "flexfoil.configuration",
    ///   "version": 1,
    ///   "ref_chord": null,
    ///   "ref_point": null,
    ///   "elements": [
    ///     {
    ///       "name": "main",
    ///       "body": { "name": "naca2412", "coordinates": [ [1.0, 0.0], ... ] },
    ///       "placement": { "pivot": [0.0, 0.0], "rotation_deg": 0.0,
    ///                      "translation": [0.0, 0.0], "scale": 1.0 },
    ///       "paneling": { "curv_param": 1.0, "te_le_ratio": 0.15,
    ///                     "te_spacing_ratio": 0.667 }
    ///     }
    ///   ]
    /// }
    /// ```
    ///
    /// `ref_chord` and `ref_point` are written as `null` when unset, so the
    /// document records that the D4 defaults apply rather than freezing the
    /// value they currently resolve to.
    ///
    /// # Errors
    /// [`GeometryError::InvalidParameter`] if any coordinate, placement or
    /// paneling value is not finite; JSON has no representation for those.
    ///
    /// # Example
    /// ```
    /// use rustfoil_core::body::Body;
    /// use rustfoil_core::configuration::Configuration;
    /// use rustfoil_core::placement::Placement;
    /// use rustfoil_core::point::point;
    ///
    /// let points = vec![
    ///     point(1.0, 0.0),
    ///     point(0.5, -0.05),
    ///     point(0.0, 0.0),
    ///     point(0.5, 0.05),
    ///     point(1.0, 0.0),
    /// ];
    /// let mut config = Configuration::single(Body::from_points("flap", &points).unwrap());
    /// config.elements[0].placement = Placement::rotation_about(point(0.7, 0.0), -30.0);
    ///
    /// let json = config.to_json().unwrap();
    /// let read_back = Configuration::from_json(&json).unwrap();
    ///
    /// assert_eq!(read_back.elements[0].placement, config.elements[0].placement);
    /// ```
    pub fn to_json(&self) -> Result<String, GeometryError> {
        let mut out = String::new();
        out.push_str("{\n  \"format\": ");
        out.push_str(&json_string(JSON_FORMAT_TAG));
        out.push_str(&format!(",\n  \"version\": {JSON_FORMAT_VERSION},\n"));

        out.push_str("  \"ref_chord\": ");
        match self.ref_chord {
            Some(chord) => out.push_str(&json_number(chord, "ref_chord")?),
            None => out.push_str("null"),
        }
        out.push_str(",\n  \"ref_point\": ");
        match self.ref_point {
            Some(p) => out.push_str(&json_pair(p.x, p.y, "ref_point")?),
            None => out.push_str("null"),
        }

        out.push_str(",\n  \"elements\": [");
        for (index, element) in self.elements.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            out.push_str("\n    {\n      \"name\": ");
            out.push_str(&json_string(&element.name));
            out.push_str(",\n      \"body\": {\n        \"name\": ");
            out.push_str(&json_string(&element.body.name));
            out.push_str(",\n        \"coordinates\": [");
            for (node_index, node) in body_contour(&element.body).iter().enumerate() {
                if node_index > 0 {
                    out.push(',');
                }
                out.push_str("\n          ");
                out.push_str(&json_pair(node.x, node.y, "coordinate")?);
            }
            out.push_str("\n        ]\n      },\n      \"placement\": ");
            out.push_str(&placement_to_json(&element.placement)?);
            out.push_str(",\n      \"paneling\": ");
            out.push_str(&paneling_to_json(&element.paneling)?);
            out.push_str("\n    }");
        }
        out.push_str("\n  ]\n}\n");

        Ok(out)
    }

    /// Read a configuration from the JSON written by
    /// [`to_json`](Configuration::to_json).
    ///
    /// `placement` and `paneling` may be omitted or `null` on an element, in
    /// which case the identity placement and [`PanelingParams::default`] apply;
    /// likewise the individual fields inside them. `format` may be omitted, so
    /// a hand-written document is not obliged to carry the tag, but a `format`
    /// that is present and different is reported rather than guessed at.
    ///
    /// # Errors
    /// [`GeometryError::ParseError`] for text that is not JSON, for a document
    /// whose shape does not match the format, and for a `format` or `version`
    /// this build cannot read. Structural problems found while walking a
    /// document that did parse report `line: 0`, since by then the position in
    /// the text is gone. Element coordinates go through
    /// [`Body::from_points`], so its errors surface unchanged.
    pub fn from_json(text: &str) -> Result<Self, GeometryError> {
        let root = JsonParser::new(text).parse_document()?;

        if let Some(tag) = root.get("format") {
            match tag.as_str() {
                Some(JSON_FORMAT_TAG) => {}
                _ => return semantic("\"format\" is not a flexfoil configuration document"),
            }
        }
        if let Some(version) = root.get("version") {
            match version.as_number() {
                Some(v) if v >= 1.0 && v <= f64::from(JSON_FORMAT_VERSION) => {}
                _ => return semantic("\"version\" is not a configuration format this build reads"),
            }
        }

        let elements_json = match root.get("elements") {
            Some(value) => match value.as_array() {
                Some(items) => items,
                None => return semantic("\"elements\" must be an array"),
            },
            None => return semantic("the document has no \"elements\" array"),
        };

        let mut elements = Vec::with_capacity(elements_json.len());
        for item in elements_json {
            elements.push(element_from_json(item)?);
        }

        let mut config = Self::new(elements);
        config.ref_chord = optional_number(&root, "ref_chord", "\"ref_chord\" must be a number")?;
        config.ref_point =
            optional_point(&root, "ref_point", "\"ref_point\" must be an [x, y] pair")?;
        Ok(config)
    }
}

// ---------------------------------------------------------------------------
// Element role names for an imported multi-element file
// ---------------------------------------------------------------------------

/// The x of an element's chord midpoint, in configuration coordinates.
fn chord_mid_x(element: &Element) -> Option<f64> {
    element
        .chord_endpoints()
        .map(|(le, te)| 0.5 * (le.x + te.x))
}

/// Name the elements of an imported multi-element configuration by their
/// geometry — see [`Configuration::from_dat_blocks`].
fn assign_default_roles(config: &mut Configuration) {
    let main = match config.main_element_index() {
        Some(index) => index,
        None => return,
    };
    let reference = match chord_mid_x(&config.elements[main]) {
        Some(x) => x,
        None => return,
    };

    let mut forward: Vec<usize> = Vec::new();
    let mut aft: Vec<usize> = Vec::new();
    for (index, element) in config.elements.iter().enumerate() {
        if index == main {
            continue;
        }
        match chord_mid_x(element) {
            Some(x) if x < reference => forward.push(index),
            _ => aft.push(index),
        }
    }

    config.elements[main].name = ROLE_MAIN.to_string();
    apply_role(&mut config.elements, &forward, ROLE_SLAT);
    apply_role(&mut config.elements, &aft, ROLE_FLAP);
}

/// Give every element in `indices` the same role, numbered in configuration
/// order when there is more than one of them.
fn apply_role(elements: &mut [Element], indices: &[usize], role: &str) {
    for (ordinal, &index) in indices.iter().enumerate() {
        elements[index].name = if indices.len() == 1 {
            role.to_string()
        } else {
            format!("{role} {}", ordinal + 1)
        };
    }
}

// ---------------------------------------------------------------------------
// `.dat` writing
// ---------------------------------------------------------------------------

/// True if a line would be read as a coordinate pair rather than as a name.
///
/// The same test the `.dat` readers apply to a file's first line.
fn reads_as_coordinates(line: &str) -> bool {
    let parts: Vec<&str> = line.split_whitespace().collect();
    parts.len() >= 2 && parts[0].parse::<f64>().is_ok() && parts[1].parse::<f64>().is_ok()
}

/// The name line for a `.dat` file — see [`Configuration::to_dat`].
fn dat_name_line(name: &str) -> String {
    let single_line: String = name
        .chars()
        .map(|c| if c == '\n' || c == '\r' { ' ' } else { c })
        .collect();

    if single_line.trim().is_empty() || reads_as_coordinates(&single_line) {
        format!("# {single_line}")
    } else {
        single_line
    }
}

/// One `.dat` coordinate line, at full `f64` precision.
fn dat_coordinate_line(p: Point) -> Result<String, GeometryError> {
    for value in [p.x, p.y] {
        if !value.is_finite() {
            return Err(GeometryError::InvalidParameter {
                name: "dat_coordinate",
                value,
            });
        }
    }
    if p.x >= DAT_SEPARATOR_SENTINEL && p.y >= DAT_SEPARATOR_SENTINEL {
        return Err(GeometryError::InvalidParameter {
            name: "dat_coordinate_reads_as_separator",
            value: p.x,
        });
    }
    Ok(format!(" {:?} {:?}", p.x, p.y))
}

// ---------------------------------------------------------------------------
// JSON writing
// ---------------------------------------------------------------------------

/// A finite `f64` as JSON: the shortest decimal that reads back identically.
fn json_number(value: f64, field: &'static str) -> Result<String, GeometryError> {
    if !value.is_finite() {
        return Err(GeometryError::InvalidParameter { name: field, value });
    }
    Ok(format!("{value:?}"))
}

/// A point or vector as a JSON `[x, y]` pair.
fn json_pair(x: f64, y: f64, field: &'static str) -> Result<String, GeometryError> {
    Ok(format!(
        "[{}, {}]",
        json_number(x, field)?,
        json_number(y, field)?
    ))
}

/// A JSON string literal, with the escapes the format requires.
fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for c in value.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// A [`Placement`] as a JSON object, on one line.
fn placement_to_json(placement: &Placement) -> Result<String, GeometryError> {
    Ok(format!(
        "{{ \"pivot\": {}, \"rotation_deg\": {}, \"translation\": {}, \"scale\": {} }}",
        json_pair(placement.pivot.x, placement.pivot.y, "placement_pivot")?,
        json_number(placement.rotation_deg, "placement_rotation_deg")?,
        json_pair(
            placement.translation.x,
            placement.translation.y,
            "placement_translation"
        )?,
        json_number(placement.scale, "placement_scale")?,
    ))
}

/// [`PanelingParams`] as a JSON object, on one line.
fn paneling_to_json(paneling: &PanelingParams) -> Result<String, GeometryError> {
    Ok(format!(
        "{{ \"curv_param\": {}, \"te_le_ratio\": {}, \"te_spacing_ratio\": {} }}",
        json_number(paneling.curv_param, "paneling_curv_param")?,
        json_number(paneling.te_le_ratio, "paneling_te_le_ratio")?,
        json_number(paneling.te_spacing_ratio, "paneling_te_spacing_ratio")?,
    ))
}

// ---------------------------------------------------------------------------
// JSON reading: text to values
// ---------------------------------------------------------------------------

/// Deepest nesting [`JsonParser`] will descend into.
///
/// The parser is recursive, so a document nested more deeply than any real
/// configuration is reported rather than allowed to exhaust the stack.
const MAX_JSON_DEPTH: usize = 32;

/// A parsed JSON value — only what the configuration format needs.
#[derive(Debug)]
enum Json {
    Null,
    /// `true` or `false`. The configuration format has no boolean field, so the
    /// value is not kept: the variant exists so that a boolean where a number
    /// or a string belongs is reported as the wrong type rather than as
    /// unrecognised syntax.
    Bool,
    Number(f64),
    Str(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    /// The value of an object member, or `None` for a missing member or a
    /// value that is not an object.
    fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(members) => members
                .iter()
                .find(|(name, _)| name == key)
                .map(|(_, value)| value),
            _ => None,
        }
    }

    fn as_number(&self) -> Option<f64> {
        match self {
            Json::Number(value) => Some(*value),
            _ => None,
        }
    }

    fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(value) => Some(value),
            _ => None,
        }
    }

    fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Array(items) => Some(items),
            _ => None,
        }
    }

    fn is_null(&self) -> bool {
        matches!(self, Json::Null)
    }

    fn is_object(&self) -> bool {
        matches!(self, Json::Object(_))
    }
}

/// A recursive-descent JSON reader over one document.
struct JsonParser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl<'a> JsonParser<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            src: text.as_bytes(),
            pos: 0,
        }
    }

    /// 1-based line number of the current position, for error reporting.
    fn line(&self) -> usize {
        let upto = self.pos.min(self.src.len());
        1 + self.src[..upto].iter().filter(|&&b| b == b'\n').count()
    }

    fn err<T>(&self, message: &'static str) -> Result<T, GeometryError> {
        Err(GeometryError::ParseError {
            line: self.line(),
            message,
        })
    }

    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    fn skip_whitespace(&mut self) {
        while let Some(b) = self.peek() {
            if matches!(b, b' ' | b'\t' | b'\n' | b'\r') {
                self.pos += 1;
            } else {
                break;
            }
        }
    }

    fn expect(&mut self, byte: u8, message: &'static str) -> Result<(), GeometryError> {
        if self.peek() == Some(byte) {
            self.pos += 1;
            Ok(())
        } else {
            self.err(message)
        }
    }

    fn literal(&mut self, word: &[u8]) -> Result<(), GeometryError> {
        if self.src[self.pos..].starts_with(word) {
            self.pos += word.len();
            Ok(())
        } else {
            self.err("unrecognised JSON literal")
        }
    }

    /// Parse one complete document, rejecting anything after it.
    fn parse_document(&mut self) -> Result<Json, GeometryError> {
        let value = self.parse_value(0)?;
        self.skip_whitespace();
        if self.pos != self.src.len() {
            return self.err("trailing characters after the JSON document");
        }
        if !value.is_object() {
            return self.err("the document is not a JSON object");
        }
        Ok(value)
    }

    fn parse_value(&mut self, depth: usize) -> Result<Json, GeometryError> {
        if depth > MAX_JSON_DEPTH {
            return self.err("the document is nested too deeply");
        }
        self.skip_whitespace();
        match self.peek() {
            Some(b'{') => self.parse_object(depth),
            Some(b'[') => self.parse_array(depth),
            Some(b'"') => Ok(Json::Str(self.parse_string()?)),
            Some(b't') => {
                self.literal(b"true")?;
                Ok(Json::Bool)
            }
            Some(b'f') => {
                self.literal(b"false")?;
                Ok(Json::Bool)
            }
            Some(b'n') => {
                self.literal(b"null")?;
                Ok(Json::Null)
            }
            Some(b'-') | Some(b'0'..=b'9') => self.parse_number(),
            Some(_) => self.err("unexpected character where a JSON value was expected"),
            None => self.err("unexpected end of input"),
        }
    }

    fn parse_object(&mut self, depth: usize) -> Result<Json, GeometryError> {
        self.expect(b'{', "expected '{'")?;
        let mut members = Vec::new();
        self.skip_whitespace();
        if self.peek() == Some(b'}') {
            self.pos += 1;
            return Ok(Json::Object(members));
        }
        loop {
            self.skip_whitespace();
            let key = self.parse_string()?;
            self.skip_whitespace();
            self.expect(b':', "expected ':' after an object key")?;
            let value = self.parse_value(depth + 1)?;
            members.push((key, value));
            self.skip_whitespace();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(Json::Object(members));
                }
                _ => return self.err("expected ',' or '}' in an object"),
            }
        }
    }

    fn parse_array(&mut self, depth: usize) -> Result<Json, GeometryError> {
        self.expect(b'[', "expected '['")?;
        let mut items = Vec::new();
        self.skip_whitespace();
        if self.peek() == Some(b']') {
            self.pos += 1;
            return Ok(Json::Array(items));
        }
        loop {
            items.push(self.parse_value(depth + 1)?);
            self.skip_whitespace();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b']') => {
                    self.pos += 1;
                    return Ok(Json::Array(items));
                }
                _ => return self.err("expected ',' or ']' in an array"),
            }
        }
    }

    fn parse_string(&mut self) -> Result<String, GeometryError> {
        self.expect(b'"', "expected '\"' at the start of a string")?;
        // Bytes are copied whole, so multi-byte characters need no special
        // handling; `\u` escapes are re-encoded as UTF-8 into the same buffer.
        let mut bytes: Vec<u8> = Vec::new();
        loop {
            let byte = match self.peek() {
                Some(byte) => byte,
                None => return self.err("unterminated string"),
            };
            self.pos += 1;
            match byte {
                b'"' => {
                    return match String::from_utf8(bytes) {
                        Ok(text) => Ok(text),
                        Err(_) => self.err("a string contains invalid UTF-8"),
                    }
                }
                b'\\' => {
                    let escape = match self.peek() {
                        Some(escape) => escape,
                        None => return self.err("unterminated escape sequence"),
                    };
                    self.pos += 1;
                    match escape {
                        b'"' => bytes.push(b'"'),
                        b'\\' => bytes.push(b'\\'),
                        b'/' => bytes.push(b'/'),
                        b'b' => bytes.push(0x08),
                        b'f' => bytes.push(0x0c),
                        b'n' => bytes.push(b'\n'),
                        b'r' => bytes.push(b'\r'),
                        b't' => bytes.push(b'\t'),
                        b'u' => {
                            let mut buffer = [0u8; 4];
                            let c = self.parse_unicode_escape()?;
                            bytes.extend_from_slice(c.encode_utf8(&mut buffer).as_bytes());
                        }
                        _ => return self.err("unrecognised escape sequence in a string"),
                    }
                }
                byte if byte < 0x20 => return self.err("unescaped control character in a string"),
                byte => bytes.push(byte),
            }
        }
    }

    fn parse_hex4(&mut self) -> Result<u32, GeometryError> {
        if self.pos + 4 > self.src.len() {
            return self.err("truncated \\u escape");
        }
        let mut value = 0u32;
        for _ in 0..4 {
            let digit = match self.src[self.pos] {
                b @ b'0'..=b'9' => u32::from(b - b'0'),
                b @ b'a'..=b'f' => u32::from(b - b'a') + 10,
                b @ b'A'..=b'F' => u32::from(b - b'A') + 10,
                _ => return self.err("invalid hex digit in a \\u escape"),
            };
            value = value * 16 + digit;
            self.pos += 1;
        }
        Ok(value)
    }

    fn parse_unicode_escape(&mut self) -> Result<char, GeometryError> {
        let first = self.parse_hex4()?;

        // A code point above the BMP arrives as a UTF-16 surrogate pair.
        if (0xD800..0xDC00).contains(&first) {
            if self.peek() != Some(b'\\') {
                return self.err("unpaired UTF-16 surrogate in a \\u escape");
            }
            self.pos += 1;
            if self.peek() != Some(b'u') {
                return self.err("unpaired UTF-16 surrogate in a \\u escape");
            }
            self.pos += 1;
            let second = self.parse_hex4()?;
            if !(0xDC00..0xE000).contains(&second) {
                return self.err("unpaired UTF-16 surrogate in a \\u escape");
            }
            let combined = 0x10000 + ((first - 0xD800) << 10) + (second - 0xDC00);
            return match char::from_u32(combined) {
                Some(c) => Ok(c),
                None => self.err("invalid code point in a \\u escape"),
            };
        }

        match char::from_u32(first) {
            Some(c) => Ok(c),
            None => self.err("invalid code point in a \\u escape"),
        }
    }

    fn parse_number(&mut self) -> Result<Json, GeometryError> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        while matches!(self.peek(), Some(b'0'..=b'9')) {
            self.pos += 1;
        }
        if self.peek() == Some(b'.') {
            self.pos += 1;
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.pos += 1;
            }
        }
        if matches!(self.peek(), Some(b'e') | Some(b'E')) {
            self.pos += 1;
            if matches!(self.peek(), Some(b'+') | Some(b'-')) {
                self.pos += 1;
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.pos += 1;
            }
        }

        // ASCII by construction, so the slice is always valid UTF-8.
        let text = core::str::from_utf8(&self.src[start..self.pos]).unwrap_or("");
        match text.parse::<f64>() {
            Ok(value) if value.is_finite() => Ok(Json::Number(value)),
            Ok(_) => self.err("a number in the document is not finite"),
            Err(_) => self.err("could not read a number"),
        }
    }
}

// ---------------------------------------------------------------------------
// JSON reading: values to a configuration
// ---------------------------------------------------------------------------

/// A problem with a document that parsed as JSON but does not match the
/// configuration format. `line: 0` because the position in the text is gone by
/// the time the shape is walked.
fn semantic<T>(message: &'static str) -> Result<T, GeometryError> {
    Err(GeometryError::ParseError { line: 0, message })
}

/// A required number.
fn number_of(value: &Json, message: &'static str) -> Result<f64, GeometryError> {
    match value.as_number() {
        Some(number) => Ok(number),
        None => semantic(message),
    }
}

/// A required `[x, y]` pair.
fn point_of(value: &Json, message: &'static str) -> Result<Point, GeometryError> {
    match value.as_array() {
        Some([x, y]) => Ok(point(number_of(x, message)?, number_of(y, message)?)),
        _ => semantic(message),
    }
}

/// An object member that may be absent or `null`, read as a number.
fn optional_number(
    object: &Json,
    key: &str,
    message: &'static str,
) -> Result<Option<f64>, GeometryError> {
    match object.get(key) {
        None => Ok(None),
        Some(value) if value.is_null() => Ok(None),
        Some(value) => Ok(Some(number_of(value, message)?)),
    }
}

/// An object member that may be absent or `null`, read as an `[x, y]` pair.
fn optional_point(
    object: &Json,
    key: &str,
    message: &'static str,
) -> Result<Option<Point>, GeometryError> {
    match object.get(key) {
        None => Ok(None),
        Some(value) if value.is_null() => Ok(None),
        Some(value) => Ok(Some(point_of(value, message)?)),
    }
}

/// One element of the `elements` array.
fn element_from_json(value: &Json) -> Result<Element, GeometryError> {
    if !value.is_object() {
        return semantic("every entry of \"elements\" must be an object");
    }

    let body_json = match value.get("body") {
        Some(body) if body.is_object() => body,
        _ => return semantic("an element has no \"body\" object"),
    };
    let body_name = match body_json.get("name") {
        None => "",
        Some(name) => match name.as_str() {
            Some(name) => name,
            None => return semantic("a body \"name\" must be a string"),
        },
    };
    let coordinates_json = match body_json.get("coordinates") {
        Some(coordinates) => match coordinates.as_array() {
            Some(items) => items,
            None => return semantic("a body's \"coordinates\" must be an array"),
        },
        None => return semantic("a body has no \"coordinates\" array"),
    };

    let mut coordinates = Vec::with_capacity(coordinates_json.len());
    for item in coordinates_json {
        coordinates.push(point_of(
            item,
            "a body coordinate must be an [x, y] pair of numbers",
        )?);
    }
    let body = Body::from_points(body_name, &coordinates)?;

    let placement = match value.get("placement") {
        None => Placement::identity(),
        Some(placement) if placement.is_null() => Placement::identity(),
        Some(placement) => placement_from_json(placement)?,
    };
    let paneling = match value.get("paneling") {
        None => PanelingParams::default(),
        Some(paneling) if paneling.is_null() => PanelingParams::default(),
        Some(paneling) => paneling_from_json(paneling)?,
    };
    let name = match value.get("name") {
        None => body.name.clone(),
        Some(name) => match name.as_str() {
            Some(name) => name.to_string(),
            None => return semantic("an element \"name\" must be a string"),
        },
    };

    Ok(Element::new(body, placement, paneling, &name))
}

/// A `placement` object; every field is optional and defaults to the identity.
fn placement_from_json(value: &Json) -> Result<Placement, GeometryError> {
    if !value.is_object() {
        return semantic("\"placement\" must be an object");
    }
    let mut placement = Placement::identity();
    if let Some(pivot) = optional_point(
        value,
        "pivot",
        "a placement \"pivot\" must be an [x, y] pair",
    )? {
        placement.pivot = pivot;
    }
    if let Some(rotation) = optional_number(
        value,
        "rotation_deg",
        "a placement \"rotation_deg\" must be a number",
    )? {
        placement.rotation_deg = rotation;
    }
    if let Some(translation) = optional_point(
        value,
        "translation",
        "a placement \"translation\" must be an [x, y] pair",
    )? {
        placement.translation = vec2(translation.x, translation.y);
    }
    if let Some(scale) = optional_number(value, "scale", "a placement \"scale\" must be a number")?
    {
        placement.scale = scale;
    }
    Ok(placement)
}

/// A `paneling` object; every field is optional and defaults to
/// [`PanelingParams::default`].
fn paneling_from_json(value: &Json) -> Result<PanelingParams, GeometryError> {
    if !value.is_object() {
        return semantic("\"paneling\" must be an object");
    }
    let mut paneling = PanelingParams::default();
    if let Some(curv_param) = optional_number(
        value,
        "curv_param",
        "a paneling \"curv_param\" must be a number",
    )? {
        paneling.curv_param = curv_param;
    }
    if let Some(te_le_ratio) = optional_number(
        value,
        "te_le_ratio",
        "a paneling \"te_le_ratio\" must be a number",
    )? {
        paneling.te_le_ratio = te_le_ratio;
    }
    if let Some(te_spacing_ratio) = optional_number(
        value,
        "te_spacing_ratio",
        "a paneling \"te_spacing_ratio\" must be a number",
    )? {
        paneling.te_spacing_ratio = te_spacing_ratio;
    }
    Ok(paneling)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    /// Split `.dat` text into a name and element coordinate blocks — the job
    /// `load_airfoil_blocks` in `rustfoil-cli` does, reduced to what the
    /// fixtures here need so that a core test does not depend on the CLI crate.
    ///
    /// Blank lines, comment lines and `999.0 999.0` sentinels are boundaries,
    /// and empty runs between consecutive boundaries are dropped. The block
    /// counts asserted below are the ones the CLI parser's own tests pin for
    /// the same files, which is what ties the two together.
    fn split_dat_blocks(text: &str) -> (String, Vec<Vec<Point>>) {
        let mut lines = text.lines().peekable();
        let first = lines.peek().copied().unwrap_or("").trim().to_string();
        let name = if reads_as_coordinates(&first) {
            String::new()
        } else {
            lines.next();
            first
        };

        let mut blocks: Vec<Vec<Point>> = Vec::new();
        let mut current: Vec<Point> = Vec::new();
        for raw_line in lines {
            let line = raw_line.trim();
            let coordinates: Option<(f64, f64)> = {
                let parts: Vec<&str> = line.split_whitespace().collect();
                match (parts.first(), parts.get(1)) {
                    (Some(x), Some(y)) => match (x.parse::<f64>(), y.parse::<f64>()) {
                        (Ok(x), Ok(y)) => Some((x, y)),
                        _ => None,
                    },
                    _ => None,
                }
            };

            let is_separator = match coordinates {
                Some((x, y)) => x >= DAT_SEPARATOR_SENTINEL && y >= DAT_SEPARATOR_SENTINEL,
                None => true,
            };
            if is_separator {
                if !current.is_empty() {
                    blocks.push(std::mem::take(&mut current));
                }
                continue;
            }
            let (x, y) = coordinates.expect("a non-separator line has coordinates");
            current.push(point(x, y));
        }
        if !current.is_empty() {
            blocks.push(current);
        }

        (name, blocks)
    }

    fn testdata(file: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../testdata")
            .join(file)
    }

    fn read_config(file: &str) -> (String, Configuration) {
        let text = std::fs::read_to_string(testdata(file)).expect("read fixture");
        let (name, blocks) = split_dat_blocks(&text);
        let config = Configuration::from_dat_blocks(&name, &blocks).expect("build configuration");
        (name, config)
    }

    /// Every element's coordinates, in configuration order.
    fn all_coordinates(config: &Configuration) -> Vec<Vec<Point>> {
        config
            .iter()
            .map(|element| body_contour(&element.body))
            .collect()
    }

    /// A closed diamond of the given chord, leading edge at `x_le`.
    fn diamond(name: &str, x_le: f64, chord: f64) -> Body {
        let points = vec![
            point(x_le + chord, 0.0),
            point(x_le + 0.5 * chord, -0.05 * chord),
            point(x_le, 0.0),
            point(x_le + 0.5 * chord, 0.05 * chord),
            point(x_le + chord, 0.0),
        ];
        Body::from_points(name, &points).unwrap()
    }

    // -- body_contour ----------------------------------------------------

    #[test]
    fn body_contour_reproduces_the_points_the_body_was_built_from() {
        // Blunt trailing edge: the contour is open, so no point is duplicated.
        let points = vec![
            point(1.0, -0.006),
            point(0.5, -0.05),
            point(0.0, 0.0),
            point(0.5, 0.05),
            point(1.0, 0.006),
        ];
        let body = Body::from_points("blunt", &points).unwrap();
        assert_eq!(body_contour(&body), points);
        assert!(!body.is_closed());

        // Sharp trailing edge: the duplicated closing point comes back too.
        let closed = diamond("sharp", 0.0, 1.0);
        let contour = body_contour(&closed);
        assert_eq!(contour.len(), closed.n_panels() + 1);
        assert_eq!(contour.first(), contour.last());

        // Round-tripping the contour rebuilds an identical body.
        let rebuilt = Body::from_points(&closed.name, &contour).unwrap();
        assert_eq!(rebuilt.n_panels(), closed.n_panels());
        assert_eq!(rebuilt.le_index(), closed.le_index());
        assert_eq!(rebuilt.is_closed(), closed.is_closed());
    }

    // -- .dat round-trips ------------------------------------------------

    /// The common case: one element, real coordinates, unchanged by a
    /// write/read cycle.
    #[test]
    fn single_element_dat_round_trip_is_exact() {
        let (name, config) = read_config("naca0012.dat");
        assert_eq!(name, "NACA 0012");
        assert_eq!(config.len(), 1);

        let written = config.to_dat(&name).unwrap();
        let (name_back, blocks_back) = split_dat_blocks(&written);
        let round_tripped = Configuration::from_dat_blocks(&name_back, &blocks_back).unwrap();

        assert_eq!(name_back, name);
        assert_eq!(all_coordinates(&round_tripped), all_coordinates(&config));
    }

    /// A single element imported from a `.dat` file must be indistinguishable
    /// from `Configuration::single`, which is what the single-element path
    /// builds today.
    #[test]
    fn one_block_matches_configuration_single() {
        let points: Vec<Point> = body_contour(&diamond("naca0012", 0.0, 1.0));
        let from_blocks =
            Configuration::from_dat_blocks("naca0012", std::slice::from_ref(&points)).unwrap();
        let direct = Configuration::single(Body::from_points("naca0012", &points).unwrap());

        assert_eq!(from_blocks.len(), 1);
        assert_eq!(from_blocks.elements[0].name, direct.elements[0].name);
        assert_eq!(
            from_blocks.elements[0].body.name,
            direct.elements[0].body.name
        );
        assert_eq!(
            from_blocks.elements[0].placement,
            direct.elements[0].placement
        );
        assert_eq!(from_blocks.ref_chord, direct.ref_chord);
    }

    /// The real thing: a McDonnell Douglas 30P-30N slat/main/flap geometry,
    /// through `.dat` -> Configuration -> `.dat` -> Configuration with every
    /// coordinate unchanged.
    #[test]
    fn real_three_element_dat_round_trip_is_exact() {
        let (name, config) = read_config("mda_30p_30n_trimmed.dat");
        assert_eq!(config.len(), 3, "fixture is a three-element geometry");

        let original = all_coordinates(&config);
        let sizes: Vec<usize> = original.iter().map(Vec::len).collect();

        let written = config.to_dat(&name).unwrap();
        let (name_back, blocks_back) = split_dat_blocks(&written);
        let round_tripped = Configuration::from_dat_blocks(&name_back, &blocks_back).unwrap();

        // Element structure survives: same count, same node counts, same order.
        assert_eq!(round_tripped.len(), 3);
        assert_eq!(
            all_coordinates(&round_tripped)
                .iter()
                .map(Vec::len)
                .collect::<Vec<_>>(),
            sizes
        );

        // Every coordinate is bit-identical, not merely close.
        assert_eq!(all_coordinates(&round_tripped), original);

        // Roles are re-derived the same way, so the second pass agrees.
        let roles: Vec<&str> = round_tripped.iter().map(|e| e.name.as_str()).collect();
        assert_eq!(roles, ["slat", "main", "flap"]);
    }

    #[test]
    fn dat_separates_elements_with_the_999_sentinel() {
        let (name, config) = read_config("mda_30p_30n_trimmed.dat");
        let written = config.to_dat(&name).unwrap();

        let separators = written
            .lines()
            .filter(|line| line.trim() == DAT_SEPARATOR_LINE.trim())
            .count();
        assert_eq!(separators, config.len() - 1);

        // The sentinel must never appear as a coordinate.
        let (_, blocks) = split_dat_blocks(&written);
        assert!(blocks
            .iter()
            .flatten()
            .all(|p| p.x < DAT_SEPARATOR_SENTINEL && p.y < DAT_SEPARATOR_SENTINEL));
    }

    #[test]
    fn dat_writes_the_geometry_a_placement_puts_on_the_page() {
        let mut config = Configuration::new(vec![
            Element::from_body(diamond("main", 0.0, 1.0)),
            Element::from_body(diamond("flap", 0.0, 0.3)),
        ]);
        config.elements[1].placement = Placement::from_translation(0.95, -0.04);

        let written = config.to_dat("placed").unwrap();
        let (_, blocks) = split_dat_blocks(&written);
        let read_back = Configuration::from_dat_blocks("placed", &blocks).unwrap();

        // The `.dat` format has nowhere to put a placement, so the placement is
        // baked into the coordinates and comes back as the identity.
        assert!(read_back.elements[1].placement.is_identity());

        // The geometry itself is unchanged: the placed flap nodes are what the
        // file carries.
        let placed: Vec<Point> = body_contour(&config.elements[1].body)
            .into_iter()
            .map(|p| config.elements[1].placement.apply(p))
            .collect();
        assert_eq!(body_contour(&read_back.elements[1].body), placed);
    }

    #[test]
    fn dat_name_line_cannot_be_mistaken_for_geometry() {
        let config = Configuration::single(diamond("main", 0.0, 1.0));

        // An ordinary name is written verbatim, as a `.dat` header always was.
        let written = config.to_dat("NACA 0012").unwrap();
        assert_eq!(written.lines().next().unwrap(), "NACA 0012");
        let (name, _) = split_dat_blocks(&written);
        assert_eq!(name, "NACA 0012");

        // A name that would read as a coordinate pair, or an empty one, is
        // commented out so that no coordinate is invented or lost.
        for awkward in ["1.0 0.0", "", "   "] {
            let written = config.to_dat(awkward).unwrap();
            let (_, blocks) = split_dat_blocks(&written);
            assert_eq!(
                blocks.iter().map(Vec::len).sum::<usize>(),
                body_contour(&config.elements[0].body).len(),
                "name {awkward:?} changed the coordinate count"
            );
        }
    }

    #[test]
    fn dat_reports_a_coordinate_it_cannot_write() {
        let mut config = Configuration::single(diamond("main", 0.0, 1.0));

        config.elements[0].placement = Placement::from_translation(f64::NAN, 0.0);
        assert!(matches!(
            config.to_dat("broken"),
            Err(GeometryError::InvalidParameter {
                name: "dat_coordinate",
                ..
            })
        ));

        // A point that would read back as an element separator is refused
        // rather than written and silently lost on the next import.
        config.elements[0].placement = Placement::from_translation(1000.0, 1000.0);
        assert!(matches!(
            config.to_dat("far away"),
            Err(GeometryError::InvalidParameter {
                name: "dat_coordinate_reads_as_separator",
                ..
            })
        ));
    }

    // -- element roles ---------------------------------------------------

    #[test]
    fn roles_come_from_the_geometry() {
        // Deliberately not in slat/main/flap file order.
        let blocks = vec![
            body_contour(&diamond("a", 0.9, 0.3)),
            body_contour(&diamond("b", -0.05, 0.15)),
            body_contour(&diamond("c", 0.0, 1.0)),
        ];
        let config = Configuration::from_dat_blocks("mixed", &blocks).unwrap();

        let roles: Vec<&str> = config.iter().map(|e| e.name.as_str()).collect();
        assert_eq!(roles, ["flap", "slat", "main"]);

        // Body names record the source file and which block each came from.
        assert_eq!(config.elements[0].body.name, "mixed element 1");
        assert_eq!(config.elements[2].body.name, "mixed element 3");
    }

    #[test]
    fn several_elements_sharing_a_role_are_numbered_in_file_order() {
        let blocks = vec![
            body_contour(&diamond("main", 0.0, 1.0)),
            body_contour(&diamond("vane", 0.9, 0.2)),
            body_contour(&diamond("flap", 1.05, 0.3)),
        ];
        let config = Configuration::from_dat_blocks("double slotted", &blocks).unwrap();

        let roles: Vec<&str> = config.iter().map(|e| e.name.as_str()).collect();
        assert_eq!(roles, ["main", "flap 1", "flap 2"]);
    }

    // -- malformed .dat input --------------------------------------------

    #[test]
    fn empty_and_malformed_dat_input_is_reported() {
        // No blocks at all.
        assert!(matches!(
            Configuration::from_dat_blocks("nothing", &[]),
            Err(GeometryError::ParseError { .. })
        ));

        // A block too short to be a contour.
        let blocks = vec![vec![point(1.0, 0.0), point(0.0, 0.0)]];
        assert!(matches!(
            Configuration::from_dat_blocks("short", &blocks),
            Err(GeometryError::InsufficientPoints {
                required: 3,
                provided: 2
            })
        ));

        // A repeated point, which would give a zero-length panel.
        let blocks = vec![vec![
            point(1.0, 0.0),
            point(0.5, 0.0),
            point(0.5, 0.0),
            point(0.0, 0.0),
        ]];
        assert!(matches!(
            Configuration::from_dat_blocks("degenerate", &blocks),
            Err(GeometryError::DegeneratePanel { .. })
        ));

        // A file that is only a header parses to no blocks, and is reported
        // rather than producing an empty configuration.
        let (name, blocks) = split_dat_blocks("JUST A HEADER\n");
        assert_eq!(name, "JUST A HEADER");
        assert!(Configuration::from_dat_blocks(&name, &blocks).is_err());
    }

    // -- JSON round-trips ------------------------------------------------

    /// A configuration whose placements, paneling and reference quantities are
    /// all non-default, so nothing can pass by luck.
    fn awkward_configuration() -> Configuration {
        let mut slat = Element::from_body(diamond("slat", -0.05, 0.15));
        slat.placement = Placement {
            pivot: point(0.02, 0.01),
            rotation_deg: 12.5,
            translation: vec2(-0.13, 0.031),
            scale: 0.97,
        };
        slat.paneling = PanelingParams {
            curv_param: 1.3,
            te_le_ratio: 0.22,
            te_spacing_ratio: 0.51,
        };
        slat.name = "slat".to_string();

        let main = Element::from_body(diamond("main", 0.0, 1.0));

        let mut flap = Element::from_body(diamond("flap", 0.0, 0.3));
        flap.placement = Placement {
            pivot: point(0.0, 0.0),
            rotation_deg: -35.75,
            translation: vec2(0.9231, -0.0187),
            scale: 1.0,
        };
        flap.paneling = PanelingParams::uniform();
        flap.name = "flap".to_string();

        Configuration::new(vec![slat, main, flap])
            .with_ref_chord(1.2345)
            .with_ref_point(point(0.31, -0.0125))
    }

    fn assert_same_configuration(read: &Configuration, original: &Configuration) {
        assert_eq!(read.len(), original.len());
        assert_eq!(read.ref_chord, original.ref_chord);
        assert_eq!(read.ref_point, original.ref_point);

        for (got, want) in read.iter().zip(original.iter()) {
            assert_eq!(got.name, want.name);
            assert_eq!(got.body.name, want.body.name);
            assert_eq!(body_contour(&got.body), body_contour(&want.body));
            // `Placement` compares exactly, which is what is wanted here.
            assert_eq!(got.placement, want.placement);
            // `PanelingParams` has no `PartialEq`; compare its fields.
            assert_eq!(got.paneling.curv_param, want.paneling.curv_param);
            assert_eq!(got.paneling.te_le_ratio, want.paneling.te_le_ratio);
            assert_eq!(
                got.paneling.te_spacing_ratio,
                want.paneling.te_spacing_ratio
            );
        }
    }

    #[test]
    fn json_round_trip_carries_placements_paneling_and_references() {
        let config = awkward_configuration();
        let json = config.to_json().unwrap();
        let read_back = Configuration::from_json(&json).unwrap();

        assert_same_configuration(&read_back, &config);

        // Non-identity placements really are non-identity on both sides.
        assert!(!read_back.elements[0].placement.is_identity());
        assert!(!read_back.elements[2].placement.is_identity());
        assert!(read_back.elements[1].placement.is_identity());

        // And a second cycle produces byte-identical text.
        assert_eq!(read_back.to_json().unwrap(), json);
    }

    #[test]
    fn json_round_trip_carries_the_real_three_element_geometry() {
        let (_, config) = read_config("mda_30p_30n_trimmed.dat");
        let json = config.to_json().unwrap();
        let read_back = Configuration::from_json(&json).unwrap();

        assert_same_configuration(&read_back, &config);
        assert_eq!(all_coordinates(&read_back), all_coordinates(&config));
    }

    #[test]
    fn json_round_trip_carries_a_single_element_configuration() {
        let (name, config) = read_config("naca0012.dat");
        let read_back = Configuration::from_json(&config.to_json().unwrap()).unwrap();

        assert_same_configuration(&read_back, &config);
        assert_eq!(read_back.elements[0].body.name, name);
    }

    #[test]
    fn json_records_unset_reference_quantities_as_null() {
        let config = Configuration::single(diamond("main", 0.0, 1.0));
        let json = config.to_json().unwrap();

        assert!(json.contains("\"ref_chord\": null"));
        assert!(json.contains("\"ref_point\": null"));

        let read_back = Configuration::from_json(&json).unwrap();
        assert_eq!(read_back.ref_chord, None);
        assert_eq!(read_back.ref_point, None);
        // The D4 defaults still resolve from the geometry.
        assert_eq!(read_back.resolved_ref_chord(), config.resolved_ref_chord());
    }

    #[test]
    fn json_defaults_a_missing_placement_and_paneling() {
        let json = r#"{
            "elements": [
              { "body": { "name": "diamond", "coordinates": [
                  [1.0, 0.0], [0.5, -0.05], [0.0, 0.0], [0.5, 0.05], [1.0, 0.0]
              ] } }
            ]
        }"#;
        let config = Configuration::from_json(json).unwrap();

        assert_eq!(config.len(), 1);
        assert!(config.elements[0].placement.is_identity());
        assert_eq!(
            config.elements[0].paneling.te_le_ratio,
            PanelingParams::default().te_le_ratio
        );
        // A missing element name falls back to the body name.
        assert_eq!(config.elements[0].name, "diamond");
    }

    #[test]
    fn json_keeps_a_name_that_needs_escaping() {
        let mut config = Configuration::single(diamond("quote\" and \\ and \ttab", 0.0, 1.0));
        config.elements[0].name = "line\nbreak \u{1f6a7} \u{7f}".to_string();

        let read_back = Configuration::from_json(&config.to_json().unwrap()).unwrap();
        assert_eq!(read_back.elements[0].name, config.elements[0].name);
        assert_eq!(
            read_back.elements[0].body.name,
            config.elements[0].body.name
        );
    }

    #[test]
    fn json_reads_escapes_it_does_not_itself_write() {
        let json = r#"{
            "elements": [
              { "name": "sl\u0061t \uD83D\uDEA7", "body": { "name": "a\/b",
                "coordinates": [[1.0, 0.0], [0.5, -0.05], [0.0, 0.0], [0.5, 0.05], [1.0, 0.0]] } }
            ]
        }"#;
        let config = Configuration::from_json(json).unwrap();

        assert_eq!(config.elements[0].name, "slat \u{1f6a7}");
        assert_eq!(config.elements[0].body.name, "a/b");
    }

    #[test]
    fn json_numbers_survive_at_full_precision() {
        // A coordinate with no short decimal form, which a fixed-precision
        // writer would round.
        let awkward = 0.1_f64 + 0.2_f64;
        let points = vec![
            point(1.0, 0.0),
            point(0.5, -awkward),
            point(0.0, 1.0 / 3.0),
            point(0.5, awkward),
            point(1.0, 0.0),
        ];
        let mut config = Configuration::single(Body::from_points("awkward", &points).unwrap());
        config.elements[0].placement.rotation_deg = -1.0 / 7.0;

        let read_back = Configuration::from_json(&config.to_json().unwrap()).unwrap();
        assert_eq!(body_contour(&read_back.elements[0].body), points);
        assert_eq!(
            read_back.elements[0].placement.rotation_deg,
            config.elements[0].placement.rotation_deg
        );
    }

    // -- malformed JSON input --------------------------------------------

    #[test]
    fn malformed_json_is_reported_not_panicked() {
        let cases = [
            "",
            "   ",
            "{",
            "}",
            "[]",
            "null",
            "42",
            "\"just a string\"",
            "{\"elements\": 3}",
            "{\"elements\": [3]}",
            "{\"elements\": [{}]}",
            // No coordinates.
            "{\"elements\": [{\"body\": {\"name\": \"a\"}}]}",
            // Coordinates that are not pairs of numbers.
            "{\"elements\": [{\"body\": {\"coordinates\": [[1.0], [0.0, 0.0], [1.0, 1.0]]}}]}",
            "{\"elements\": [{\"body\": {\"coordinates\": [[\"a\", 0.0]]}}]}",
            // Unterminated structures and strings.
            "{\"elements\": [",
            "{\"elements\": [{\"body\": {\"name\": \"unterminated",
            "{\"format\": \"flexfoil.configuration\", \"version\": 1,}",
            // Trailing content after a complete document.
            "{\"elements\": []} trailing",
            // Bad numbers and escapes.
            "{\"ref_chord\": 1.0e, \"elements\": []}",
            "{\"ref_chord\": -, \"elements\": []}",
            "{\"elements\": [], \"note\": \"\\q\"}",
            "{\"elements\": [], \"note\": \"\\u00\"}",
            "{\"elements\": [], \"note\": \"\\uD83D only\"}",
            // Wrong types for the reference quantities.
            "{\"elements\": [], \"ref_chord\": \"one\"}",
            "{\"elements\": [], \"ref_point\": [0.0]}",
            // A format or version this build does not read.
            "{\"format\": \"something.else\", \"elements\": []}",
            "{\"format\": \"flexfoil.configuration\", \"version\": 99, \"elements\": []}",
            // Placement and paneling of the wrong shape.
            "{\"elements\": [{\"body\": {\"coordinates\": [[1.0, 0.0], [0.0, 0.0], [1.0, 1.0]]}, \"placement\": 3}]}",
            "{\"elements\": [{\"body\": {\"coordinates\": [[1.0, 0.0], [0.0, 0.0], [1.0, 1.0]]}, \"paneling\": {\"curv_param\": \"lots\"}}]}",
        ];

        for case in cases {
            let result = Configuration::from_json(case);
            assert!(
                result.is_err(),
                "expected an error for {case:?}, got a configuration"
            );
        }
    }

    #[test]
    fn deeply_nested_json_is_reported_rather_than_exhausting_the_stack() {
        let depth = MAX_JSON_DEPTH * 4;
        let nested = format!(
            "{{\"elements\": {}{}}}",
            "[".repeat(depth),
            "]".repeat(depth)
        );
        assert!(matches!(
            Configuration::from_json(&nested),
            Err(GeometryError::ParseError { .. })
        ));
    }

    #[test]
    fn an_empty_configuration_round_trips_as_json() {
        let config = Configuration::new(vec![]);
        let json = config.to_json().unwrap();
        let read_back = Configuration::from_json(&json).unwrap();

        assert!(read_back.is_empty());
        assert_eq!(read_back.to_json().unwrap(), json);
    }

    #[test]
    fn json_reports_a_value_it_cannot_write() {
        let mut config = Configuration::single(diamond("main", 0.0, 1.0));
        config.elements[0].placement.scale = f64::INFINITY;
        assert!(matches!(
            config.to_json(),
            Err(GeometryError::InvalidParameter {
                name: "placement_scale",
                ..
            })
        ));

        let mut config = Configuration::single(diamond("main", 0.0, 1.0));
        config.elements[0].paneling.curv_param = f64::NAN;
        assert!(matches!(
            config.to_json(),
            Err(GeometryError::InvalidParameter {
                name: "paneling_curv_param",
                ..
            })
        ));

        let config = Configuration::single(diamond("main", 0.0, 1.0)).with_ref_chord(f64::NAN);
        assert!(matches!(
            config.to_json(),
            Err(GeometryError::InvalidParameter {
                name: "ref_chord",
                ..
            })
        ));
    }
}
