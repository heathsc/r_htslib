use super::{bam1_t, bcf1_t, tbx_t, get_cstr, hts_err, kstring_t, sam_hdr_t, SamHeader};

use c2rust_bitfields::BitfieldStruct;
use libc::{c_char, c_int, c_short, c_uchar, c_uint};
use std::ffi::{c_void, CString};
use std::io;
use std::marker::PhantomData;
use std::ptr::{null_mut, NonNull};
use std::ops::{Deref, DerefMut};

pub type HtsPos = i64;

pub const HTS_IDX_NOCOOR: HtsPos = -2;
pub const HTS_IDX_START: HtsPos = -3;
pub const HTS_IDX_REST: HtsPos = -4;
pub const HTS_IDX_NONE: HtsPos = -5;

pub const FT_UNKN: u32 = 0;
pub const FT_GZ: u32 = 1;
pub const FT_VCF: u32 = 2;
pub const FT_VCF_GZ: u32 = FT_GZ | FT_VCF;
pub const FT_BCF: u32 = 4;
pub const FT_BCF_GZ: u32 = FT_GZ | FT_BCF;
pub const FT_STDIN: u32 = 8;

#[repr(C)]
pub struct hts_opt {
    _unused: [u8; 0],
}

#[repr(C)]
pub struct hts_filter_t {
    _unused: [u8; 0],
}

#[repr(C)]
#[derive(BitfieldStruct)]
pub struct htsFile {
    #[bitfield(name = "is_bin", ty = "c_uchar", bits = "0..=0")]
    #[bitfield(name = "is_write", ty = "c_uchar", bits = "1..=1")]
    #[bitfield(name = "is_be", ty = "c_uchar", bits = "2..=2")]
    #[bitfield(name = "is_cram", ty = "c_uchar", bits = "3..=3")]
    #[bitfield(name = "is_bgzf", ty = "c_uchar", bits = "4..=4")]
    #[bitfield(name = "dummy", ty = "u32", bits = "5..=31")]
    bfield: [u8; 4],
    lineno: i64,
    line: kstring_t,
    fn_: *mut c_char,
    fn_aux: *mut c_char,
    fp: *mut c_void,
    state: *mut c_void,
    format: htsFormat,
    idx: *mut hts_idx_t,
    fnidx: *const c_char,
    bam_header: *mut sam_hdr_t,
    filter: *mut hts_filter_t,
}

impl htsFile {
    pub fn format(&self) -> HtsFormat {
        HtsFormat {
            inner: NonNull::new(unsafe { hts_get_format(self) }).expect("hts_get_format returned NULL"),
            phantom: PhantomData,
            owned: false,
        }
    }
    /*
    pub fn sam_index_load(&mut self) -> io::Result<HtsIndex> {
        match NonNull::new(unsafe {
            sam_index_load(self.inner_mut(), get_cstr(&self.name).as_ptr())
        }) {
            None => Err(hts_err(format!(
                "Couldn't load SAM index for file {}",
                self.name
            ))),
            Some(idx) => Ok(HtsIndex {
                inner: idx,
                phantom: PhantomData,
            }),
        }
    } */
}

#[repr(C)]
pub struct hts_idx_t {
    _unused: [u8; 0],
}

#[repr(C)]
pub struct hts_reglist_t {
    _unused: [u8; 0],
}

#[repr(C)]
pub struct hts_pair64_max {
    _unused: [u8; 0],
}

#[repr(C)]
pub struct HtsItrBins {
    n: c_int,
    m: c_int,
    a: *mut c_int,
}

#[repr(C)]
#[derive(BitfieldStruct)]
pub struct hts_itr_t {
    #[bitfield(name = "read_rest", ty = "c_uchar", bits = "0..=0")]
    #[bitfield(name = "finished", ty = "c_uchar", bits = "1..=1")]
    #[bitfield(name = "is_cram", ty = "c_uchar", bits = "2..=2")]
    #[bitfield(name = "nocoor", ty = "c_uchar", bits = "3..=3")]
    #[bitfield(name = "multi", ty = "c_uchar", bits = "4..=4")]
    #[bitfield(name = "dummy", ty = "u32", bits = "5..=31")]
    bfield: [u8; 4],
    tid: c_int,
    n_off: c_int,
    i: c_int,
    n_reg: c_int,
    beg: HtsPos,
    end: HtsPos,
    reg_list: *mut hts_reglist_t,
    curr_tid: c_int,
    curr_reg: c_int,
    curr_intv: c_int,
    curr_beg: HtsPos,
    curr_end: HtsPos,
    curr_off: u64,
    nocoor_off: u64,
    off: *mut hts_pair64_max,
    readrec: HtsReadrecFunc,
    seek: HtsSeekFunc,
    tell: HtsTellFunc,
    bins: HtsItrBins,
}

pub type HtsReadrecFunc = unsafe extern fn(fp: *mut BGZF, data: *mut c_void, r: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos);
pub type HtsSeekFunc = unsafe extern fn(fp: *mut c_void, offset: i64, where_: c_int);
pub type HtsTellFunc = unsafe extern fn(fp: *mut c_void);

#[repr(C)]
pub struct BGZF {
    _unused: [u8; 0],
}
#[repr(C)]
pub enum htsFormatCategory {
    UnknownCategory,
    SequenceData,
    VariantData,
    IndexFile,
    RegionList,
}

impl Default for htsFormatCategory {
    fn default() -> Self {
        Self::UnknownCategory
    }
}

#[repr(C)]
#[derive(PartialEq)]
pub enum htsExactFormat {
    UnknownFormat,
    BinaryFormat,
    TextFormat,
    Sam,
    Bam,
    Bai,
    Cram,
    Crai,
    Vcf,
    Bcf,
    Csi,
    Gzi,
    Tbi,
    Bed,
    HtsGet,
    EmptyFormat,
    FastaFormat,
    FastqFormat,
    FaiFormat,
    FqiFormat,
    HtsCrypt4GH,
}

impl Default for htsExactFormat {
    fn default() -> Self {
        Self::UnknownFormat
    }
}

#[repr(C)]
#[derive(PartialEq)]
pub enum htsCompression {
    NoCompression,
    Gzip,
    Bgzf,
    Custom,
    Bzip2Compression,
}

impl Default for htsCompression {
    fn default() -> Self {
        Self::NoCompression
    }
}

#[repr(C)]
#[derive(Default)]
pub struct htsFormatVersion {
    major: c_short,
    minor: c_short,
}

#[repr(C)]
pub struct htsFormat {
    category: htsFormatCategory,
    format: htsExactFormat,
    version: htsFormatVersion,
    compression: htsCompression,
    compression_level: c_short,
    specific: *mut hts_opt,
    phantom: PhantomData<hts_opt>,
}

impl Default for htsFormat {
    fn default() -> Self {
        Self {
            category: htsFormatCategory::default(),
            format: htsExactFormat::default(),
            version: htsFormatVersion::default(),
            compression: htsCompression::default(),
            compression_level: 0,
            specific: std::ptr::null_mut::<hts_opt>(),
            phantom: PhantomData,
        }
    }
}

impl Drop for htsFormat {
    fn drop(&mut self) {
        if !self.specific.is_null() {
            unsafe { hts_opt_free(self.specific) }
        }
    }
}

#[link(name = "hts")]
extern "C" {
    fn hts_open(fn_: *const c_char, mode: *const c_char) -> *mut htsFile;
    fn hts_open_format(fn_: *const c_char, mode: *const c_char, fmt: *const htsFormat) -> *mut htsFile;
    fn hts_close(fp_: *mut htsFile) -> c_int;
    fn hts_set_threads(fp_: *mut htsFile, t_: c_int) -> c_int;
    fn hts_get_format(fp_: *const htsFile) -> *mut htsFormat;
    fn sam_index_load(fp_: *mut htsFile, name: *const c_char) -> *mut hts_idx_t;
    fn hts_set_fai_filename(fp_: *mut htsFile, fn_aux: *const c_char) -> c_int;
    fn hts_itr_destroy(iter: *mut hts_itr_t);
    fn sam_itr_queryi(idx: *const hts_idx_t, tid: c_int, start: HtsPos, end: HtsPos) -> *mut hts_itr_t;
    fn sam_itr_regarray(idx: *const hts_idx_t, hdr: *mut sam_hdr_t, regarray: *const *const c_char, count: c_uint) -> *mut hts_itr_t;
    fn hts_itr_multi_next(fp: *mut htsFile, itr: *mut hts_itr_t, r: *mut c_void) -> c_int;
    fn hts_itr_next(fp: *mut BGZF, itr: *mut hts_itr_t, r: *mut c_void, data: *mut c_void) -> c_int;
    fn hts_parse_format(format: *mut htsFormat, str: *const c_char) -> c_int;
    fn hts_opt_free(opts: *mut hts_opt);
    fn hts_opt_add(opts: *mut *mut hts_opt, c_arg: *const c_char) -> c_int;
}

pub struct HtsFile {
    inner: NonNull<htsFile>,
    name: String,
    phantom: PhantomData<htsFile>,
}

unsafe impl Send for HtsFile {}

impl Deref for HtsFile {
    type Target = htsFile;
    #[inline]
    fn deref(&self) -> &htsFile { unsafe{self.inner.as_ref()} }
}

impl DerefMut for HtsFile {
    #[inline]
    fn deref_mut(&mut self) -> &mut htsFile {unsafe{ self.inner.as_mut() }}
}

impl AsRef<htsFile> for HtsFile {
    #[inline]
    fn as_ref(&self) -> &htsFile { self}
}

impl AsMut<htsFile> for HtsFile {
    #[inline]
    fn as_mut(&mut self) -> &mut htsFile { self}
}

impl Drop for HtsFile {
    fn drop(&mut self) {
        unsafe { hts_close(self.as_mut()) };
    }
}

impl HtsFile {
    pub fn new<S: AsRef<str>>(name: S, mode: &str) -> io::Result<Self> {
        let name = name.as_ref();
        let hfile = unsafe { hts_open(get_cstr(name).as_ptr(), get_cstr(mode).as_ptr()) };
        Self::from_hts_file(hfile, name)
    }

    pub fn open_format<S: AsRef<str>>(
        name: S,
        mode: &[c_char],
        fmt: &HtsFormat,
    ) -> io::Result<Self> {
        let name = name.as_ref();
        let hfile =
            unsafe { hts_open_format(get_cstr(name).as_ptr(), mode.as_ptr(), fmt.inner.as_ptr()) };
        Self::from_hts_file(hfile, name)
    }

    fn from_hts_file(hfile: *mut htsFile, name: &str) -> io::Result<Self> {
        match NonNull::new(hfile) {
            None => Err(hts_err(format!("Couldn't open file {}", name))),
            Some(fptr) => Ok(HtsFile {
                inner: fptr,
                name: name.to_owned(),
                phantom: PhantomData,
            }),
        }
    }

    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    pub fn set_threads(&mut self, t: usize) -> io::Result<()> {
        let ret = unsafe { hts_set_threads(self.as_mut(), t as c_int) };
        if ret != 0 {
            Err(hts_err(format!(
                "Failed to set additional threads to file {}",
                self.name
            )))
        } else {
            Ok(())
        }
    }

    pub fn sam_index_load(&mut self) -> io::Result<HtsIndex> {
        match NonNull::new(unsafe {
            sam_index_load(self.as_mut(), get_cstr(&self.name).as_ptr())
        }) {
            None => Err(hts_err(format!(
                "Couldn't load SAM index for file {}",
                self.name
            ))),
            Some(idx) => Ok(HtsIndex {
                inner: idx,
                phantom: PhantomData,
            }),
        }
    }

    pub fn set_fai_filename<S: AsRef<str>>(&mut self, name: S) -> io::Result<()> {
        let name = name.as_ref();
        let ret = unsafe { hts_set_fai_filename(self.as_mut(), get_cstr(name).as_ptr()) };
        if ret != 0 {
            Err(hts_err(format!(
                "Failed to attach reference index {} to file {}",
                name, self.name
            )))
        } else {
            Ok(())
        }
    }
}

pub struct HtsIndex {
    inner: NonNull<hts_idx_t>,
    phantom: PhantomData<hts_idx_t>,
}

unsafe impl Send for HtsIndex {}

impl HtsIndex {
    fn inner(&self) -> &hts_idx_t {
        unsafe { self.inner.as_ref() }
    }
    pub fn sam_itr_queryi(&self, tid: isize, start: usize, end: usize) -> io::Result<HtsItr> {
        let it = NonNull::new(unsafe {
            sam_itr_queryi(
                self.inner(),
                tid as libc::c_int,
                start as HtsPos,
                end as HtsPos,
            )
        });
        if let Some(itr) = it {
            Ok(HtsItr {
                inner: itr,
                itr_type: HtsItrType::SamItr,
                phantom: PhantomData,
            })
        } else {
            Err(hts_err("Failed to obtain sam iterator".to_string()))
        }
    }
    pub fn sam_itr_regarray(&self, hdr: &mut SamHeader, regions: &[String]) -> io::Result<HtsItr> {
        let count = regions.len();
        if count == 0 {
            return self.sam_itr_queryi(HTS_IDX_START as isize, 0, 0);
        }
        // We need to do this in 2 stages: generate and array of CStrings and then an array of ptrs to the CStrings
        // otherwise the CStrings go out of scope before they are used and the ptrs point to garbage
        let carray: Vec<CString> = regions.iter().map(get_cstr).collect();
        let parray: Vec<*const c_char> = carray.iter().map(|cs| cs.as_ptr()).collect();
        let it = NonNull::new(unsafe {
            sam_itr_regarray(
                self.inner(),
                hdr.as_mut(),
                parray.as_ptr(),
                count as c_uint,
            )
        });
        if let Some(itr) = it {
            Ok(HtsItr {
                inner: itr,
                itr_type: HtsItrType::SamItr,
                phantom: PhantomData,
            })
        } else {
            Err(hts_err("Failed to obtain sam iterator".to_string()))
        }
    }
}

pub struct HtsFormat {
    inner: NonNull<htsFormat>,
    phantom: PhantomData<htsFormat>,
    owned: bool,
}

impl Drop for HtsFormat {
    fn drop(&mut self) {
        if self.owned {
            let fmt = unsafe { Box::from_raw(self.inner.as_ptr()) };
            drop(fmt)
        }
    }
}

impl Default for HtsFormat {
    fn default() -> Self {
        let fmt = Box::new(htsFormat::default());
        let fmt: *mut htsFormat = Box::into_raw(fmt);
        Self {
            inner: NonNull::new(fmt).unwrap(),
            phantom: PhantomData,
            owned: true,
        }
    }
}

impl HtsFormat {
    pub fn parse_format<S: AsRef<str>>(s: S) -> io::Result<Self> {
        let s = s.as_ref();
        let fmt = Box::new(htsFormat::default());
        let fmt: *mut htsFormat = Box::into_raw(fmt);
        if unsafe { hts_parse_format(fmt, get_cstr(s).as_ptr()) } == 0 {
            Ok(Self {
                inner: NonNull::new(fmt).unwrap(),
                phantom: PhantomData,
                owned: true,
            })
        } else {
            Err(hts_err(format!("Error parsing '{}' as hts format", s)))
        }
    }

    pub fn opt_add<S: AsRef<str>>(&mut self, s: S) -> io::Result<()> {
        let s = s.as_ref();
        if unsafe { hts_opt_add(&mut self.inner.as_mut().specific, get_cstr(s).as_ptr()) } != 0 {
            Err(hts_err(format!("hts_add_opt(): Illegal option {}", s)))
        } else {
            Ok(())
        }
    }

    pub fn is_compressed(&self) -> bool {
        let format = unsafe { &self.inner.as_ref() };
        format.compression == htsCompression::Bgzf || format.format == htsExactFormat::Cram
    }
}

pub enum HtsItrRec {
    BamRec(bam1_t),
    KStr(kstring_t),
    BcfRec(bcf1_t),
}

pub enum HtsItrType<'a> {
    SamItr,
    TbxItr(&'a mut tbx_t),
    VcfItr,
    None,
}

impl <'a> HtsItrType<'a> {
    pub fn take(&mut self) -> HtsItrType<'a> {
        std::mem::replace(self, Self::None)
    }
}

pub struct HtsItr<'a> {
    inner: NonNull<hts_itr_t>,
    itr_type: HtsItrType<'a>,
    phantom: PhantomData<hts_itr_t>,
}

impl <'a>Deref for HtsItr<'a> {
    type Target = hts_itr_t;
    #[inline]
    fn deref(&self) -> &hts_itr_t { unsafe{self.inner.as_ref()} }
}

impl <'a>DerefMut for HtsItr<'a> {
    #[inline]
    fn deref_mut(&mut self) -> &mut hts_itr_t {unsafe{ self.inner.as_mut() }}
}

impl <'a>AsRef<hts_itr_t> for HtsItr<'a> {
    #[inline]
    fn as_ref(&self) -> &hts_itr_t { self}
}

impl <'a>AsMut<hts_itr_t> for HtsItr<'a> {
    #[inline]
    fn as_mut(&mut self) -> &mut hts_itr_t { self}
}

unsafe impl <'a> Send for HtsItr<'a> {}

impl <'a>Drop for HtsItr<'a> {
    fn drop(&mut self) {
        unsafe { hts_itr_destroy(self.inner.as_ptr()) };
    }
}

impl <'a> HtsItr<'a> {
    pub fn new(itr: *mut hts_itr_t, itr_type: HtsItrType<'a>) -> Option<Self> {
        if matches!(itr_type, HtsItrType::None) {
            None
        } else {
            NonNull::new(itr).map(|p| HtsItr { inner: p, itr_type, phantom: PhantomData })
        }
    }

    pub fn next(&mut self, fp: &mut HtsFile, rec: &mut HtsItrRec) -> io::Result<bool> {

       let i = if matches!(self.itr_type, HtsItrType::SamItr) && self.as_ref().multi() != 0 {
            if let HtsItrRec::BamRec(b) = rec {
                unsafe { hts_itr_multi_next(fp.as_mut(), self.as_mut(), b as *mut bam1_t as *mut c_void) }
            } else {
                return Err(hts_err("Wrong record type for HtsItr::next()".to_string()))
            }
        } else {
            let bgfp = if fp.as_ref().is_bgzf() != 0 {
                fp.as_ref().fp as *mut BGZF
            } else {
                null_mut::<BGZF>()
            };
            match (&mut self.itr_type, rec) {
                (HtsItrType::SamItr, HtsItrRec::BamRec(b)) => {
                    unsafe {hts_itr_next(bgfp, self.as_mut(), b as *mut bam1_t as *mut c_void, fp.as_mut() as *mut htsFile as *mut c_void)}
                },
                (HtsItrType::VcfItr, HtsItrRec::BcfRec(b)) if !bgfp.is_null() => {
                    unsafe {hts_itr_next(bgfp, self.as_mut(), b as *mut bcf1_t as *mut c_void, null_mut())}
                },
                (HtsItrType::TbxItr(_), HtsItrRec::KStr(s)) if !bgfp.is_null() => {
                    let mut t = self.itr_type.take();
                    let tbx = if let HtsItrType::TbxItr(tx) = &mut t { tx } else { panic!("Wrong HtsItrType") };
                    let j = unsafe {hts_itr_next(bgfp, self.as_mut(), s as *mut kstring_t as *mut c_void, *tbx as *mut tbx_t as *mut c_void)};
                    self.itr_type = t;
                    j
                },
                (HtsItrType::VcfItr, _) | (HtsItrType::TbxItr(_), _) if bgfp.is_null() => return Err(hts_err("Only bgzf compressed files can be used with iterators".to_string())),
                _ => return Err(hts_err("Wrong record type for HtsItr::next()".to_string())),
            }
        };
        if i >= 0 {
            Ok(true)
        } else if i == -1 {
            Ok(false)
        } else {
            Err(hts_err("Error reading record".to_string()))
        }
    }
/*
    pub fn sam_itr_next(&mut self, fp: &mut HtsFile, mut brec: BamRec) -> SamReadResult {
        if self.itr_type != HtsItrType::SamItr {
            return SamReadResult::Error
        }
        let p = brec.as_mut();
        match unsafe {
            if (*self.inner.as_ref()).multi() != 0 {
                hts_itr_multi_next(
                    fp.inner_mut(),
                    self.inner.as_ptr(),
                    p as *mut bam1_t as *mut c_void,
                )
            } else {
                hts_itr_next(
                    if fp.inner().is_bgzf() != 0 {
                        fp.inner().fp as *mut BGZF
                    } else {
                        null_mut::<BGZF>()
                    },
                    self.inner.as_ptr(),
                    p as *mut bam1_t as *mut c_void,
                    fp.inner_mut() as *mut htsFile as *mut c_void,
                )
            }
        } {
            0..=c_int::MAX => SamReadResult::Ok(brec),
            -1 => SamReadResult::EOF,
            _ => SamReadResult::Error,
        }
    } */
}
