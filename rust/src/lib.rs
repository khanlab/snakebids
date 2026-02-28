//! Rust-accelerated internals for snakebids.
//!
//! This module exposes `parse_format_string`, a fast equivalent of the Python
//! `SnakemakeFormatter.parse()` parsing loop.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

const UNEXPECTED_CLOSE: &str = "unexpected '}' in string";
const MISSING_CLOSE: &str = "expected '}' before end of string";
const UNEXPECTED_OPEN: &str = "unexpected '{' in field name";

/// Parse a Snakemake-style format string character by character.
///
/// Returns a list of 5-tuples per parsed segment:
///   `(literal_text, field_name, new_underscore, new_current_field, update_current_field)`
///
/// - `literal_text`         – the literal portion before the next field (or at end of string)
/// - `field_name`           – `None` for non-field segments; the field name for field segments
/// - `new_underscore`       – `None` = no change; `Some(s)` = set `_underscore` to `s`
/// - `new_current_field`    – new `_current_field` value (only meaningful when `update_current_field`)
/// - `update_current_field` – whether `_current_field` should be updated
///
/// Raises `ValueError` on malformed input (same conditions as the Python implementation).
#[pyfunction]
pub fn parse_format_string(
    format_string: &str,
) -> PyResult<Vec<(String, Option<String>, Option<String>, Option<String>, bool)>> {
    let mut entries: Vec<(String, Option<String>, Option<String>, Option<String>, bool)> =
        Vec::new();

    let mut chars = format_string.chars();

    // Accumulates literal text between fields.
    let mut literal = String::new();

    loop {
        match chars.next() {
            // ---- End of string ------------------------------------------
            None => {
                if !literal.is_empty() {
                    entries.push((literal, None, None, None, false));
                }
                return Ok(entries);
            }

            // ---- Opening brace ------------------------------------------
            Some('{') => {
                match chars.next() {
                    None => {
                        // Trailing lone `{` — no closing brace
                        return Err(PyValueError::new_err(MISSING_CLOSE));
                    }
                    Some('{') => {
                        // `{{` — escaped open brace (emits up-to-and-including first `{`)
                        literal.push('{');
                        entries.push((literal, None, Some("_".to_string()), None, false));
                        literal = String::new();
                    }
                    Some('}') => {
                        // `{}` — empty field name
                        let new_underscore = literal.chars().last().map(|last| {
                            if last == '/' || last == '_' {
                                String::new()
                            } else {
                                "_".to_string()
                            }
                        });
                        entries.push((
                            literal,
                            Some(String::new()),
                            new_underscore,
                            None,
                            true,
                        ));
                        literal = String::new();
                    }
                    Some(first) => {
                        // Start of a real field.  Collect everything up to the matching `}`.
                        let mut field_content = String::new();
                        field_content.push(first);

                        let closing = loop {
                            match chars.next() {
                                None => {
                                    return Err(PyValueError::new_err(MISSING_CLOSE));
                                }
                                Some('{') => {
                                    // Nested `{` inside a field is not allowed.
                                    return Err(PyValueError::new_err(UNEXPECTED_OPEN));
                                }
                                Some('}') => break field_content,
                                Some(c) => field_content.push(c),
                            }
                        };

                        // Split field name from constraint at the first `,`.
                        let (field_name, new_current_field) =
                            match closing.find(',') {
                                Some(comma) => {
                                    let name = closing[..comma].to_string();
                                    let full = closing.clone();
                                    (name, Some(full))
                                }
                                None => (closing, None),
                            };

                        // `_underscore` update based on the last char of the literal.
                        let new_underscore: Option<String> =
                            literal.chars().last().map(|last| {
                                if last == '/' || last == '_' {
                                    String::new()
                                } else {
                                    "_".to_string()
                                }
                            });

                        entries.push((
                            literal,
                            Some(field_name),
                            new_underscore,
                            new_current_field,
                            true,
                        ));
                        literal = String::new();
                    }
                }
            }

            // ---- Closing brace ------------------------------------------
            Some('}') => {
                match chars.next() {
                    Some('}') => {
                        // `}}` — escaped close brace
                        literal.push('}');
                        entries.push((literal, None, Some("_".to_string()), None, false));
                        literal = String::new();
                    }
                    _ => {
                        // Lone `}` outside a field
                        return Err(PyValueError::new_err(UNEXPECTED_CLOSE));
                    }
                }
            }

            // ---- Ordinary character -------------------------------------
            Some(c) => {
                literal.push(c);
            }
        }
    }
}

/// Register this module as `snakebids._rust._core`.
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_format_string, m)?)?;
    Ok(())
}
