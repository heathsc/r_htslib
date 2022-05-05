use std::io;
use std::marker::PhantomData;
use std::ptr::NonNull;
use std::ops::{Deref, DerefMut};
use std::convert::{AsRef, AsMut};

pub mod lib;
pub use lib::*;

use super::{
    from_cstr, get_cstr, htsFile, hts_err, kstring_t, Hts, HtsFile, HtsPos, HtsHdr,
    HtsRead, HtsItr, BGZF, hts_itr_next, HtsFileDesc
};
use libc::{c_int, c_void};

#[derive(Debug)]
pub struct VcfHeader {
    inner: NonNull<bcf_hdr_t>,
    phantom: PhantomData<bcf_hdr_t>,
}

impl Deref for VcfHeader {
    type Target = bcf_hdr_t;
    #[inline]
    fn deref(&self) -> &bcf_hdr_t { unsafe {self.inner.as_ref()} }
}

impl DerefMut for VcfHeader {
    #[inline]
    fn deref_mut(&mut self) -> &mut bcf_hdr_t { unsafe {self.inner.as_mut() }}
}

impl AsRef<bcf_hdr_t> for VcfHeader {
    #[inline]
    fn as_ref(&self) -> &bcf_hdr_t { self}
}

impl AsMut<bcf_hdr_t> for VcfHeader {
    #[inline]
    fn as_mut(&mut self) -> &mut bcf_hdr_t { self}
}

impl VcfHeader {
    pub fn new(mode: &str) -> io::Result<VcfHeader> {
        match NonNull::new(unsafe { bcf_hdr_init(get_cstr(mode).as_ptr()) }) {
            None => Err(hts_err("Couldn't create VCF/BCF header".to_string())),
            Some(hdr) => Ok(VcfHeader { inner: hdr , phantom: PhantomData}),
        }
    }

    pub fn read(fp: &mut HtsFile) -> io::Result<VcfHeader> {
        match NonNull::new(unsafe{bcf_hdr_read(fp.as_mut())}) {
            None => Err(hts_err(format!("Failed to load header from {}", fp.name()))),
            Some(p) => {
                Ok(Self {
                    inner: p,
                    phantom: PhantomData,
                })
            },
        }
    }

    pub fn append<S: AsRef<str>>(&mut self, line: S) -> io::Result<()> {
        match unsafe { bcf_hdr_append(self.as_mut(), get_cstr(line).as_ptr()) } {
            0 => Ok(()),
            _ => Err(hts_err(
                "Error appending line to VCF/BCF header".to_string(),
            )),
        }
    }

    pub fn get_version(&self) -> &str {
        from_cstr(unsafe { bcf_hdr_get_version(self.as_ref()) })
    }

    pub fn add_sample<S: AsRef<str>>(&mut self, name: S) -> io::Result<()> {
        match unsafe { bcf_hdr_add_sample(self.as_mut(), get_cstr(name).as_ptr()) } {
            0 => Ok(()),
            _ => Err(hts_err("Error adding sample to VCF/BCF header".to_string())),
        }
    }
    pub fn write(&mut self, hout: &mut HtsFile) -> io::Result<()> {
        match unsafe { bcf_hdr_write(hout.as_mut(), self.as_mut()) } {
            0 => Ok(()),
            _ => Err(hts_err("Error writing VCF/BCF header".to_string())),
        }
    }
    pub fn id2int<S: AsRef<str>>(&self, category: usize, name: S) -> Option<usize> {
        let ix = unsafe {
            bcf_hdr_id2int(
                self.as_ref(),
                category as c_int,
                get_cstr(name.as_ref()).as_ptr(),
            )
        };
        if ix < 0 {
            None
        } else {
            Some(ix as usize)
        }
    }
    pub fn sync(&mut self) -> io::Result<()> {
        match unsafe { bcf_hdr_sync(self.as_mut()) } {
            0 => Ok(()),
            _ => Err(hts_err("Error adding sample to VCF/BCF header".to_string())),
        }
    }
}

impl Drop for VcfHeader {
    fn drop(&mut self) {
        unsafe { bcf_hdr_destroy(self.as_mut()) };
    }
}


pub struct BcfRec {
    inner: NonNull<bcf1_t>,
    phantom: PhantomData<bcf1_t>,
    line: kstring_t,
}

unsafe impl Send for BcfRec {}

impl Deref for BcfRec {
    type Target = bcf1_t;
    #[inline]
    fn deref(&self) -> &bcf1_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for BcfRec {
    #[inline]
    fn deref_mut(&mut self) -> &mut bcf1_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<bcf1_t> for BcfRec {
    #[inline]
    fn as_ref(&self) -> &bcf1_t { self}
}

impl AsMut<bcf1_t> for BcfRec {
    #[inline]
    fn as_mut(&mut self) -> &mut bcf1_t { self}
}

impl HtsRead for BcfRec {
    fn read(&mut self, hts: &mut Hts) -> io::Result<bool> {
        let (fp, hdr) = hts.hts_file_and_header();
        let hts_file = fp.as_mut();
        let res = if let Some(HtsHdr::Vcf(hd)) = &hdr {
            let i = unsafe { bcf_read(hts_file, hd.as_ref(), self.as_mut()) };
            match i {
                0..=c_int::MAX if self.errcode == 0 => Ok(true),
                -1 => Ok(false),
                _ => Err(hts_err("Error reading VCF/BCF record".to_string())),
            }
        } else {
            Err(hts_err(format!("Appropriate header not found for file {}", fp.name())))
        };
        res
    }

    fn read_itr(&mut self, hts: &mut Hts, itr: &mut HtsItr) -> io::Result<bool> {
        let (fp, hdr) = hts.hts_file_and_header();
        if let Some(HtsFileDesc::Bgzf(bgzf)) = fp.file_desc() {
            if let Some(HtsHdr::Vcf(h)) = hdr {
                match unsafe { hts_itr_next(
                    bgzf.as_ptr(),
                    itr.as_mut(),
                    self as *mut BcfRec as *mut c_void,
                    h.as_mut() as *mut bcf_hdr_t as *mut c_void,
                )} {
                    0..=c_int::MAX => Ok(true),
                    -1 => Ok(false),
                    _ => Err(hts_err("Error reading VCF/BCF record".to_string())),
                }
            } else {
                Err(hts_err(format!("Appropriate header not found for file {}", fp.name())))
            }
        } else {
            Err(hts_err(format!("File {} is not in bgzf format (required for indexing)", fp.name())))
        }
    }
}

impl BcfRec {
    pub fn new() -> io::Result<Self> {
        match NonNull::new(unsafe { bcf_init() }) {
            Some(b) => Ok(Self { inner: b, phantom: PhantomData, line: kstring_t::new() }),
            None => Err(hts_err("Failed to allocate new BcfRec".to_string())),
        }
    }

    pub fn set_n_allele(&mut self, n_all: usize) {
        self.as_mut().set_n_allele(n_all as u16)
    }
    pub fn set_n_info(&mut self, n_info: usize) {
        self.as_mut().set_n_info(n_info as u16)
    }
    pub fn set_n_sample(&mut self, n_sample: usize) {
        self.as_mut().set_n_sample(n_sample as u32)
    }
    pub fn set_n_fmt(&mut self, n_fmt: usize) {
        self.as_mut().set_n_fmt(n_fmt as u8)
    }
    pub fn set_qual(&mut self, qual: f32) {
        self.as_mut().qual = qual
    }
    pub fn write(&mut self, file: &mut HtsFile, hdr: &mut VcfHeader) -> io::Result<()> {
        if unsafe { bcf_write(file.as_mut(), hdr.as_mut(), self.as_mut()) } < 0 {
            Err(hts_err("Error writing out VCF record".to_string()))
        } else {
            Ok(())
        }
    }
    pub fn format(&mut self, hdr: &HtsHdr, s: &mut kstring_t) -> io::Result<()>{
        s.clear();
        if let HtsHdr::Vcf(h) = hdr {
            let e = self.as_mut().format(h, s);
            if e >= 0 {
                Ok(())
            } else {
                Err(hts_err(format!("Error formatting VCF record: err {}", e)))
            }
        } else {
            Err(hts_err("Wrong header type for VCF format".to_string()))
        }
    }
}

impl Drop for BcfRec {
    fn drop(&mut self) {
        unsafe { bcf_destroy(self.inner.as_ptr()) }
    }
}

