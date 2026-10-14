import React, {useState} from 'react';

export default function PolarProvenanceInspector() {
  const [metadata, setMetadata] = useState(true);
  const [matching, setMatching] = useState(true);
  return <fieldset><legend>Distribution provenance</legend>
    <label><input type="checkbox" checked={metadata} onChange={e=>setMetadata(e.target.checked)} /> Installed distribution metadata exists</label>{' '}
    <label><input type="checkbox" checked={matching} onChange={e=>setMatching(e.target.checked)} /> Metadata identifies this loaded module</label>
    <p>solver_version: {metadata && matching ? 'the matching installed release' : 'None (unknown)'}. A polar can still be generated in either case. This demonstration does not query your environment or establish native binary identity.</p>
  </fieldset>;
}
