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
/// Returns a list of 4-tuples per parsed segment:
///   `(literal_text, field_name, squelch_underscore, constraint)`
///
/// - `literal_text`       – the literal portion before the next field (or at end of string)
/// - `field_name`         – `None` for non-field segments; the field name for field segments
/// - `squelch_underscore` – `None` = no change to `_underscore` (empty literal before a field,
///                          or an end-of-string segment); `Some(true)` = set `_underscore = ""`
///                          (literal ends in `'/'` or `'_'`); `Some(false)` = set `_underscore = "_"`
///                          (all other non-empty literals, including doubled-brace segments)
/// - `constraint`         – `None` for non-field segments; `""` for a field with no constraint;
///                          `",<text>"` (starting with `,`) for a field that has a constraint;
///                          Python can derive `_current_field` as `field_name + constraint`
///
/// Raises `ValueError` on malformed input (same conditions as the Python implementation).
#[pyfunction]
pub fn parse_format_string(
    format_string: &str,
) -> PyResult<Vec<(String, Option<String>, Option<bool>, Option<String>)>> {
    let mut entries: Vec<(String, Option<String>, Option<bool>, Option<String>)> = Vec::new();

    let mut chars = format_string.chars();

    // Accumulates literal text between fields.
    let mut literal = String::new();

    loop {
        match chars.next() {
            // ---- End of string ------------------------------------------
            None => {
                if !literal.is_empty() {
                    entries.push((literal, None, None, None));
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
                        // `{{` — escaped open brace; always sets _underscore to "_"
                        literal.push('{');
                        entries.push((literal, None, Some(false), None));
                        literal = String::new();
                    }
                    Some(first) => {
                        // `{}` or `{name}` or `{name,constraint}` — a real field.
                        // Collect everything up to the matching `}`.
                        let mut field_content = String::new();
                        if first != '}' {
                            field_content.push(first);
                            loop {
                                match chars.next() {
                                    None => {
                                        return Err(PyValueError::new_err(MISSING_CLOSE));
                                    }
                                    Some('{') => {
                                        // Nested `{` inside a field is not allowed.
                                        return Err(PyValueError::new_err(UNEXPECTED_OPEN));
                                    }
                                    Some('}') => break,
                                    Some(c) => field_content.push(c),
                                }
                            }
                        }

                        // Split field name from constraint at the first `,`.
                        let (field_name, constraint) = match field_content.find(',') {
                            Some(comma) => {
                                let name = field_content[..comma].to_string();
                                let cons = field_content[comma..].to_string();
                                (name, cons)
                            }
                            None => (field_content, String::new()),
                        };

                        // `squelch_underscore`: None when literal is empty (no update needed),
                        // Some(true) when literal ends in a squelcher ('/' or '_'),
                        // Some(false) otherwise.
                        let squelch = literal.chars().last().map(|c| c == '/' || c == '_');

                        entries.push((literal, Some(field_name), squelch, Some(constraint)));
                        literal = String::new();
                    }
                }
            }

            // ---- Closing brace ------------------------------------------
            Some('}') => {
                match chars.next() {
                    Some('}') => {
                        // `}}` — escaped close brace; always sets _underscore to "_"
                        literal.push('}');
                        entries.push((literal, None, Some(false), None));
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
