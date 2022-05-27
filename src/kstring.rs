use std::ptr::NonNull;
use std::marker::PhantomData;

use super::*;
use libc::{c_char, c_float, c_int, c_void, size_t};

#[repr(C)]
#[derive(Debug)]
pub struct kstring_t {
    l: size_t,
    m: size_t,
    s: Option<NonNull<c_char>>,
    phantom: PhantomData<kstring_t>,
}

impl Default for kstring_t {
    fn default() -> Self {
        Self {
            l: 0,
            m: 0,
            s: None,
            phantom: PhantomData,
        }
    }
}

impl Drop for kstring_t {
    fn drop(&mut self) {
        if let Some(p) = self.s {
            unsafe { libc::free(p.as_ptr() as *mut c_void) }
        }
    }
}

#[link(name = "hts")]
extern "C" {
    fn bcf_enc_vchar(s: *mut kstring_t, l: c_int, a: *const c_char) -> c_int;
    fn bcf_enc_vint(s: *mut kstring_t, n: c_int, a: *const i32, wsize: c_int) -> c_int;
    fn bcf_enc_vfloat(s: *mut kstring_t, l: c_int, a: *const c_float) -> c_int;
}

impl kstring_t {
    pub fn new() -> Self { Self::default() }

    pub fn len(&self) -> size_t { self.l }

    pub fn clear(&mut self) { self.l = 0 }

    pub fn resize(&mut self, size: size_t) -> bool {
        if self.m < size {
            let size = if size > (usize::MAX >> 2) {
                size
            } else {
                size + (size >> 1)
            };
            let p = NonNull::new(if let Some(ptr) = &mut self.s {
                unsafe { libc::realloc(ptr.as_ptr() as *mut c_void, size) }
            } else {
                unsafe { libc::malloc(size) }
            }).map(|x| x.cast::<c_char>());
            if p.is_none() { return true }
            self.s = p;
            self.m = size;
        }
        false
    }

    pub fn putsn(&mut self, p: *const c_char, l: size_t) -> bool {
        let new_sz = self.l + l + 2;
        if new_sz <= self.l || self.resize(new_sz) {
            true
        } else {
            let l1 = self.l as isize;
            unsafe {
                let ptr = self.s.unwrap().as_ptr().offset(l1);
                libc::memcpy(ptr as *mut c_void, p as *const c_void, l);
                self.l += l;
                *(ptr.offset(l as isize)) = 0;
            }
            false
        }
    }

    pub fn putc(&mut self, c: c_char) -> bool {
        if self.resize(self.l + 2) {
            true
        } else {
            let l = self.l as isize;
            unsafe {
                let ptr = self.s.unwrap().as_ptr();
                *(ptr.offset(l)) = c;
                *(ptr.offset(l + 1)) = 0;
            }
            self.l += 1;
            false
        }
    }

    pub fn bcf_enc_size(&mut self, size: c_int, bcf_type: u8) -> bool {
        if size >= 15 {
            self.putc((15 << 4 | bcf_type) as c_char)
                || if size >= 128 {
                    if size >= 32768 {
                        self.putc((1 << 4 | BCF_BT_INT32) as c_char)
                            || self
                                .putsn((size as c_int).to_le_bytes().as_ptr() as *const c_char, 4)
                    } else {
                        self.putc((1 << 4 | BCF_BT_INT16) as c_char)
                            || self.putsn((size as u16).to_le_bytes().as_ptr() as *const c_char, 2)
                    }
                } else {
                    self.putc((1 << 4 | BCF_BT_INT8) as c_char) || self.putc(size as c_char)
                }
        } else {
            self.putc(((size as u8) << 4 | bcf_type) as c_char)
        }
    }
    pub fn bcf_enc_int1(&mut self, x: c_int) -> bool {
        if x == BCF_INT32_VECTOR_END {
            self.bcf_enc_size(1, BCF_BT_INT8 as u8) || self.putc(BCF_INT8_VECTOR_END as c_char)
        } else if x == BCF_INT32_MISSING {
            self.bcf_enc_size(1, BCF_BT_INT8 as u8) || self.putc(BCF_INT8_MISSING as c_char)
        } else if x <= BCF_MAX_BT_INT8 && x >= BCF_MIN_BT_INT8 {
            self.bcf_enc_size(1, BCF_BT_INT8 as u8) || self.putc(x as c_char)
        } else if x <= BCF_MAX_BT_INT16 && x >= BCF_MIN_BT_INT16 {
            self.bcf_enc_size(1, BCF_BT_INT16 as u8)
                || self.putsn((x as u16).to_le_bytes().as_ptr() as *const i8, 2)
        } else {
            self.bcf_enc_size(1, BCF_BT_INT32 as u8)
                || self.putsn(x.to_le_bytes().as_ptr() as *const i8, 4)
        }
    }
    pub fn bcf_enc_vchar(&mut self, v: &[u8]) -> bool {
        unsafe {
            bcf_enc_vchar(
                self as *mut kstring_t,
                v.len() as c_int,
                v.as_ptr() as *const c_char,
            ) < 0
        }
    }
    pub fn bcf_enc_vint(&mut self, v: &[i32]) -> bool {
        unsafe { bcf_enc_vint(self as *mut kstring_t, v.len() as c_int, v.as_ptr(), -1) < 0 }
    }
    pub fn bcf_enc_vfloat(&mut self, v: &[f32]) -> bool {
        unsafe { bcf_enc_vfloat(self as *mut kstring_t, v.len() as c_int, v.as_ptr()) < 0 }
    }
    pub fn to_str(&self) -> Option<&str> {
        self.s.map(|s| from_cstr(s.as_ptr()).trim_end())
    }
    pub fn to_cstr(&self) -> Option<&CStr> {
        self.s.map(|s| unsafe { CStr::from_ptr(s.as_ptr())})
    }
    pub fn as_ptr(&self) -> Option<NonNull<c_char>> { self.s }

    pub fn as_slice(&self, inc_zero: bool) -> Option<&[u8]> {
        self.s.map(|s| {
            let p = s.as_ptr() as *const u8;
            unsafe {std::slice::from_raw_parts(p, if inc_zero { self.l +  1 } else { self.l })}
        })
    }
    pub fn as_slice_mut(&mut self, inc_zero: bool) -> Option<&mut [u8]> {
        self.s.map(|s| {
            let p = s.as_ptr() as *mut u8;
            unsafe {std::slice::from_raw_parts_mut(p, if inc_zero { self.l +  1 } else { self.l })}
        })
    }
}
