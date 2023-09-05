#[macro_use]
extern crate anyhow;

use std::{
    ffi::{CStr, CString},
    io::{self, Error, ErrorKind},
};

use anyhow::Context;

pub mod hts;
pub use hts::*;
pub mod sam;
pub use sam::*;
pub mod vcf;
pub use vcf::*;
pub mod faidx;
pub use faidx::*;
pub mod kstring;
pub mod tbx;
pub use tbx::*;
pub mod regidx;
pub use kstring::*;

/// Make a CString from &str.  Panics if input &str contains NUL bytes
#[inline]
pub fn get_cstr<S: AsRef<str>>(s: S) -> CString {
    CString::new(s.as_ref().as_bytes()).unwrap()
}

/// Make a &str from a C string ptr. Panics if input ptr is null or if string contain non
/// UTF8 sequences
#[inline]
pub fn from_cstr<'a>(c: *const i8) -> &'a str {
    c_to_cstr(c).to_str().expect("String not UTF8")
}

/// Try to make a &str from a C string ptr.  Returns Error if ptr is null of if string contains non
/// UTF8 sequences
#[inline]
pub fn try_from_cstr<'a>(c: *const i8) -> anyhow::Result<&'a str> {
    if c.is_null() {
        Err(anyhow!("Null ptr in try_from_cstr()"))
    } else {
        unsafe {
            CStr::from_ptr(c)
                .to_str()
                .with_context(|| "Conversion error in try_from_cstr()")
        }
    }
}

/// Convert from &str to CStr.  Panics if input ptr is null
#[inline]
pub fn c_to_cstr<'a>(c: *const i8) -> &'a CStr {
    if c.is_null() {
        panic!("from_cstr() called with a NULL");
    }
    unsafe { CStr::from_ptr(c) }
}

pub fn hts_err(s: String) -> io::Error {
    Error::new(ErrorKind::Other, s)
}
