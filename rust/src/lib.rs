//! Rust-accelerated internals for snakebids.
//!
//! This module exposes `parse_format_string`, a fast equivalent of the Python
//! `SnakemakeFormatter.parse()` parsing loop.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

const UNEXPECTED_CLOSE: &str = "unexpected '}' in string";
const MISSING_CLOSE: &str = "expected '}' before end of string";
const UNEXPECTED_OPEN: &str = "unexpected '{' in field name";

/// Find `target` in `chars[start..]`, returning an absolute index.
#[inline]
fn find_char(chars: &[char], target: char, start: usize) -> Option<usize> {
    chars[start..]
        .iter()
        .position(|&c| c == target)
        .map(|p| start + p)
}

/// Find `target` in `chars[start..end]`, returning an absolute index.
#[inline]
fn find_char_in(chars: &[char], target: char, start: usize, end: usize) -> Option<usize> {
    chars[start..end]
        .iter()
        .position(|&c| c == target)
        .map(|p| start + p)
}

/// Parse a Snakemake-style format string.
///
/// Mirrors the logic of `SnakemakeFormatter.parse()` exactly, including all side
/// effects on `_underscore` and `_current_field`.
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
    let chars: Vec<char> = format_string.chars().collect();
    let length = chars.len();

    let mut entries: Vec<(String, Option<String>, Option<String>, Option<String>, bool)> =
        Vec::new();

    // Mirrors Python locals: i=-1, anchor=0, next_char='{', j=first '}'
    let mut i: Option<usize> = None; // None ≡ Python's -1
    let mut anchor: usize = 0;
    let mut next_char: char = '{';
    let mut j: usize = find_char(&chars, '}', 0).unwrap_or(length);

    while anchor < length {
        // i = format_string.find(next_char, i + 1)
        let search_start = i.map(|v| v + 1).unwrap_or(0);
        i = find_char(&chars, next_char, search_start);

        // Mirrors: if i > j: swap / elif i == -1: ...
        match i {
            Some(i_val) if i_val > j => {
                i = Some(j);
                j = i_val;
            }
            None => {
                if j == length {
                    // No braces left – emit trailing literal and return.
                    let literal: String = chars[anchor..].iter().collect();
                    entries.push((literal, None, None, None, false));
                    return Ok(entries);
                }
                i = Some(j);
                j = length;
            }
            _ => {}
        }

        let i_val = i.unwrap(); // safe: set above

        // if final character
        if i_val + 1 == length {
            if chars[i_val] == '{' {
                return Err(PyValueError::new_err(MISSING_CLOSE));
            }
            return Err(PyValueError::new_err(UNEXPECTED_CLOSE));
        }

        // if doubled brace ( {{ or }} )
        if chars[i_val + 1] == chars[i_val] {
            let i_new = i_val + 1;
            // literal up to (and including) first brace of the pair
            let literal: String = chars[anchor..i_new].iter().collect();
            entries.push((literal, None, Some("_".to_string()), None, false));
            anchor = i_new + 1;
            next_char = chars[i_new];
            i = Some(i_new);
            continue;
        }

        // if undoubled '}' outside a field
        if chars[i_val] == '}' {
            return Err(PyValueError::new_err(UNEXPECTED_CLOSE));
        }

        // ---- Opening brace: parse field --------------------------------
        let literal: String = chars[anchor..i_val].iter().collect();

        if j == length {
            return Err(PyValueError::new_err(MISSING_CLOSE));
        }

        // No nested '{' allowed between i+1 and j
        if find_char_in(&chars, '{', i_val + 1, j).is_some() {
            return Err(PyValueError::new_err(UNEXPECTED_OPEN));
        }

        // Extract field name (strip constraint after first ',')
        let (field_name, new_current_field) =
            match find_char_in(&chars, ',', i_val + 1, j) {
                Some(comma) => {
                    let name: String = chars[i_val + 1..comma].iter().collect();
                    let full: String = chars[i_val + 1..j].iter().collect();
                    (name, Some(full))
                }
                None => {
                    let name: String = chars[i_val + 1..j].iter().collect();
                    (name, None)
                }
            };

        // _underscore update: based on literal_text's last character
        let new_underscore: Option<String> = if !literal.is_empty() {
            let last = literal.chars().last().unwrap();
            if last == '/' || last == '_' {
                Some(String::new())
            } else {
                Some("_".to_string())
            }
        } else {
            None // no change when literal is empty
        };

        // Advance j to the next '}' after the current closing brace
        let j_val = j;
        j = find_char(&chars, '}', j_val + 1).unwrap_or(length);

        anchor = j_val + 1;
        next_char = '{';
        i = Some(j_val);

        entries.push((literal, Some(field_name), new_underscore, new_current_field, true));
    }

    Ok(entries)
}

/// Register this module as `snakebids._rust._core`.
#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(parse_format_string, m)?)?;
    Ok(())
}
