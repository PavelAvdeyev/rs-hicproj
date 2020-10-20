use std::{error, fmt};

#[derive(Debug, Clone)]
pub struct MatrixIndexError;

impl fmt::Display for MatrixIndexError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "The incorrect index is provided.")
    }
}

impl error::Error for MatrixIndexError {}

#[derive(Debug, Clone)]
pub struct SelectorUninitError;

impl fmt::Display for SelectorUninitError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Selector for matrix is not initialized.")
    }
}

impl error::Error for SelectorUninitError {}

#[derive(Debug, Clone)]
pub struct MatrixResolutionError;

impl fmt::Display for MatrixResolutionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Such resolution does not exist.")
    }
}

impl error::Error for MatrixResolutionError {}
