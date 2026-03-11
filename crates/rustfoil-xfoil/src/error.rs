use thiserror::Error;

#[derive(Debug, Error)]
pub enum XfoilError {
    #[error("{0}")]
    Message(String),
    #[error("geometry error: {0}")]
    Geometry(String),
    #[error("inviscid solver error: {0}")]
    Inviscid(String),
}

pub type Result<T> = std::result::Result<T, XfoilError>;

impl From<rustfoil_inviscid::InviscidError> for XfoilError {
    fn from(value: rustfoil_inviscid::InviscidError) -> Self {
        Self::Inviscid(value.to_string())
    }
}
