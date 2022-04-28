use super::{bam1_t, bcf1_t, tbx_t, get_cstr, from_cstr, hts_err, kstring_t, sam_hdr_t};

use c2rust_bitfields::BitfieldStruct;
use libc::{c_char, c_int, c_short, c_uchar, c_uint};
use std::ffi::{c_void, CString};
use std::io;
use std::marker::PhantomData;
use std::ptr::{null_mut, null, NonNull};
use std::ops::{Deref, DerefMut};

use crate::{
    tbx_index_load3, VcfHeader, SamHeader, Tbx, tbx_readrec,
    vcf_read_itr, bcf_read_itr, tbx_read_itr
};

pub type HtsPos = i64;

pub const HTS_IDX_NOCOOR: HtsPos = -2;
pub const HTS_IDX_START: HtsPos = -3;
pub const HTS_IDX_REST: HtsPos = -4;
pub const HTS_IDX_NONE: HtsPos = -5;

pub const HTS_IDX_SILENT_FAIL: c_int = 2;

pub const HTS_FMT_CSI: c_int = 0;
pub const HTS_FMT_BAI: c_int = 1;
pub const HTS_FMT_TBI: c_int = 2;
pub const HTS_FMT_CRAI: c_int = 3;
pub const HTS_FMT_FAI: c_int = 4;

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
struct hts_tpool {
    _unused: [u8; 0],
}

#[repr(C)]
pub struct prehtsFile {
    bfield: [u8; 4],
    lineno: i64,
    line: kstring_t,
    fn_: *mut c_char,
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
    fn_: NonNull<c_char>,
    fn_aux: *mut c_char,
    pub (crate) fp: *mut c_void,
    state: *mut c_void,
    format: htsFormat,
    idx: *mut hts_idx_t,
    fnidx: *const c_char,
    bam_header: *mut sam_hdr_t,
    filter: *mut hts_filter_t,
}

pub enum HtsFileDesc {
    Bgzf(NonNull<BGZF>),
    Cram(NonNull<cram_fd>),
    Hfile(NonNull<hfile>),
}

impl htsFile {
    pub fn format(&self) -> &htsFormat { &self.format }

    pub fn name(&self) -> &str {
        from_cstr(self.name_ptr())
    }

    pub fn file_desc(&mut self) -> Option<HtsFileDesc> {
        if self.is_bgzf() != 0 {
            NonNull::new(self.fp as *mut BGZF).map(|p| HtsFileDesc::Bgzf(p))
        } else if self.is_cram() != 0 {
            NonNull::new(self.fp as *mut cram_fd).map(|p| HtsFileDesc::Cram(p))
        } else {
            NonNull::new(self.fp as *mut hfile).map(|p| HtsFileDesc::Hfile(p))
        }
    }

    fn name_ptr(&self) -> *const c_char {
        self.fn_.as_ptr()
    }

    pub fn set_threads(&mut self, t: usize) -> io::Result<()> {
        let ret = unsafe { hts_set_threads(self, t as c_int) };
        if ret != 0 {
            Err(hts_err(format!("Failed to set threads for file {}", self.name())))
        } else {
            Ok(())
        }
    }

    pub fn set_fai_filename<S: AsRef<str>>(&mut self, name: S) -> io::Result<()> {
        let name = name.as_ref();
        let ret = unsafe { hts_set_fai_filename(self, get_cstr(name).as_ptr()) };
        if ret != 0 {
            Err(hts_err(format!("Failed to attach reference index {} to file {}", name, self.name())))
        } else {
            Ok(())
        }
    }
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

pub type HtsReadrecFunc = unsafe extern fn(fp: *mut BGZF, data: *mut c_void, r: *mut c_void, tid: *mut c_int, beg: *mut HtsPos, end: *mut HtsPos) -> c_int;
pub type HtsSeekFunc = unsafe extern fn(fp: *mut c_void, offset: i64, where_: c_int);
pub type HtsTellFunc = unsafe extern fn(fp: *mut c_void);

#[repr(C)]
pub struct BGZF {
    _unused: [u8; 0],
}

#[repr(C)]
pub struct cram_fd {
    _unused: [u8; 0],
}
#[repr(C)]
pub struct hfile {
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
    UnknownFormat, BinaryFormat, TextFormat,
    Sam, Bam, Bai, Cram, Crai, Vcf, Bcf, Csi, Gzi, Tbi, Bed, HtsGet,
    EmptyFormat, FastaFormat, FastqFormat, FaiFormat, FqiFormat, HtsCrypt4GH,
}

impl Default for htsExactFormat {
    fn default() -> Self { Self::UnknownFormat }
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

impl htsFormat {
    pub fn format(&self) -> &htsExactFormat { &self.format }
}

#[link(name = "hts")]
extern "C" {
    fn hts_open_format(fn_: *const c_char, mode: *const c_char, fmt: *const htsFormat) -> *mut prehtsFile;
    fn hts_close(fp_: *mut htsFile) -> c_int;
    fn hts_set_threads(fp_: *mut htsFile, t_: c_int) -> c_int;
    fn hts_set_thread_pool(fp_: *mut htsFile, p: *const HtsThreadPool) -> c_int;
    fn hts_tpool_init(n: c_int) -> *mut hts_tpool;
    fn hts_tpool_destroy(p: *mut hts_tpool);
    fn sam_index_load3(fp_: *mut htsFile, name: *const c_char, fnidx: *const c_char, flags: c_int) -> *mut hts_idx_t;
    fn hts_idx_load3(name: *const c_char, fnidx: *const c_char, fmt: c_int, flags : c_int) -> *mut hts_idx_t;
    fn hts_idx_destroy(idx: *mut hts_idx_t);
    fn hts_set_fai_filename(fp_: *mut htsFile, fn_aux: *const c_char) -> c_int;
    fn hts_itr_destroy(iter: *mut hts_itr_t);
    fn hts_itr_query(idx: *const hts_idx_t, tid: c_int, begin: HtsPos, end: HtsPos, readrec: HtsReadrecFunc) -> *mut hts_itr_t;
    fn sam_itr_queryi(idx: *const hts_idx_t, tid: c_int, start: HtsPos, end: HtsPos) -> *mut hts_itr_t;
//    fn sam_itr_regarray(idx: *const hts_idx_t, hdr: *mut sam_hdr_t, regarray: *const *const c_char, count: c_uint) -> *mut hts_itr_t;
    pub (crate) fn hts_itr_multi_next(fp: *mut htsFile, itr: *mut hts_itr_t, r: *mut c_void) -> c_int;
    pub(crate) fn hts_itr_next(fp: *mut BGZF, itr: *mut hts_itr_t, r: *mut c_void, data: *mut c_void) -> c_int;
    fn hts_parse_format(format: *mut htsFormat, str: *const c_char) -> c_int;
    fn hts_opt_free(opts: *mut hts_opt);
    fn hts_opt_add(opts: *mut *mut hts_opt, c_arg: *const c_char) -> c_int;
    pub(crate) fn bgzf_getline(fp: *mut BGZF, delim: c_int, str: *mut kstring_t) -> c_int;
}

pub struct Hts {
    pub(crate) hts_file: HtsFile,
    pub(crate) header: Option<HtsHdr>,
    pub(crate) idx: Option<HtsIndex>,
}

unsafe impl Send for Hts {}

impl Hts {
    pub fn open<S: AsRef<str>>(name: S, mode: &str) -> io::Result<Self> {
        Self::open_format_(name, mode, None)
    }

    pub fn open_format<S: AsRef<str>>(name: S, mode: &str, fmt: &HtsFormat) -> io::Result<Self> {
        Self::open_format_(name, mode, Some(fmt))
    }

    fn open_format_<S: AsRef<str>>(name: S, mode: &str, fmt: Option<&HtsFormat>) -> io::Result<Self> {
        let name = name.as_ref();
        let mut fp = HtsFile::open_format_(name, mode, fmt)?;
        let header = if fp.as_ref().is_write() == 0 {
            match fp.format().format {
                htsExactFormat::Sam | htsExactFormat::Bam | htsExactFormat::Cram => Some(HtsHdr::Sam(SamHeader::read(&mut fp)?)),
                htsExactFormat::Bcf | htsExactFormat::Vcf => Some(HtsHdr::Vcf(VcfHeader::read(&mut fp)?)),
                _ => if let Ok(tbx) = Tbx::load(name) { Some(HtsHdr::Tbx(tbx)) } else { None }
            }
        } else { None };
        Ok(Self {
            hts_file: fp,
            header,
            idx: None,
        })
    }

    pub fn header(&self) -> Option<&HtsHdr> { self.header.as_ref() }

    pub fn header_mut(&mut self) -> Option<&mut HtsHdr> { self.header.as_mut() }

    pub fn hts_file(&mut self) -> &mut HtsFile { &mut self.hts_file }

    pub fn hts_file_and_header(&mut self) -> (&mut HtsFile, Option<&mut HtsHdr>) {
        (&mut self.hts_file, self.header.as_mut())
    }

    pub fn index_load(&mut self) -> io::Result<()> {
        self.index_load3(None, 0)
    }

    pub fn index_load2(&mut self, fnidx: Option<&str>) -> io::Result<()> {
        self.index_load3(fnidx, 0)
    }

    pub fn name(&self) -> &str { self.hts_file.name() }

    pub fn has_index(&self) -> bool {
        self.idx.is_some() || matches!(self.header, Some(HtsHdr::Tbx(_)))
    }

    pub fn index_load3(&mut self, fnidx: Option<&str>, flags: usize) -> io::Result<()> {
        let flags = (flags as c_int) | HTS_IDX_SILENT_FAIL;
        let fname = self.hts_file.name_ptr();
        let fnidx = fnidx.map(get_cstr);
        let have_fnidx = fnidx.is_some();
        let fnidx_ptr = match fnidx {
            Some(s) => s.as_ptr(),
            None => null(),
        };

        let mk_hts_idx = |p: *mut hts_idx_t| {
            NonNull::new(p).map(|q| HtsIndex { inner: q, phantom: PhantomData})
        };

        let idx = match self.hts_file.format().format {
            htsExactFormat::Bam | htsExactFormat::Cram => {
                self.idx = mk_hts_idx(unsafe { sam_index_load3(self.hts_file.as_mut(), fname, fnidx_ptr, flags) });
                self.idx.is_some()
            },
            htsExactFormat::Sam => {
                self.idx = mk_hts_idx(unsafe { sam_index_load3(self.hts_file.as_mut(), fname, fnidx_ptr, flags) });
                self.idx.is_some() || {
                    if let Some(q) = NonNull::new(unsafe { tbx_index_load3(fname, fnidx_ptr, flags) }) {
                        let tbx = Tbx::new(q);
                        self.header = Some(HtsHdr::Tbx(tbx));
                        true
                    } else { false }
                }
            },
            htsExactFormat::Bcf => {
                self.idx = mk_hts_idx(unsafe { hts_idx_load3(fname, fnidx_ptr, HTS_FMT_CSI, flags)});
                self.idx.is_some()
            },
            htsExactFormat::Vcf => {
                self.idx = mk_hts_idx(unsafe { hts_idx_load3(fname, fnidx_ptr, HTS_FMT_CSI, flags)});
                self.idx.is_some() || {
                    if let Some(q) = NonNull::new(unsafe { tbx_index_load3(fname, fnidx_ptr, flags) }) {
                        let tbx = Tbx::new(q);
                        self.header = Some(HtsHdr::Tbx(tbx));
                        true
                    } else { false }
                }
            },
            _ => {
                // Retry to load index for Tabix file if we have been supplied with an index name
                if have_fnidx {
                    if let Some(q) = NonNull::new(unsafe { tbx_index_load3(fname, fnidx_ptr, flags) }) {
                        let tbx = Tbx::new(q);
                        self.header = Some(HtsHdr::Tbx(tbx));
                    }
                }
                self.header.is_some()
            }
        };
        if idx { Ok(()) } else { Err(hts_err(format!("Failed to load index for file {}", self.name()))) }
    }

    pub fn index(&self) -> Option<&hts_idx_t> {
        if self.has_index() {
            self.idx.as_ref().map_or_else(
                || if let Some(HtsHdr::Tbx(tbx)) = self.header.as_ref() {
                    Some(tbx.idx())
                } else { None },
                |p| Some(p.as_ref()))
        } else { None }
    }

    pub fn rec_type(&self) -> Option<HtsRecType> {
        self.header().map(|h| match h {
            HtsHdr::Sam(_) => HtsRecType::Sam,
            HtsHdr::Vcf(_) => HtsRecType::Vcf,
            HtsHdr::Tbx(_) => HtsRecType::Tbx,
        })
    }

    pub fn itr_queryi(&mut self, tid: usize, begin: HtsPos, end: HtsPos) -> io::Result<HtsItr> {

        // Load index if not already present
        if !self.has_index() { self.index_load()? }

        // Extract index
        let idx = self.index().unwrap();

        // And header
        let hdr = self.header.as_ref().unwrap();

        if let Some(itr) = NonNull::new(match hdr {
            HtsHdr::Sam(_) => unsafe { sam_itr_queryi(idx, tid as c_int, begin, end) },
            HtsHdr::Tbx(_) => unsafe { hts_itr_query(idx, tid as c_int, begin, end, tbx_read_itr) },
            HtsHdr::Vcf(_) if matches!(self.hts_file.format().format, htsExactFormat::Bcf) => unsafe { hts_itr_query(idx, tid as c_int, begin, end, bcf_read_itr) },
            HtsHdr::Vcf(_) => unsafe { hts_itr_query(idx, tid as c_int, begin, end, vcf_read_itr) },
        }) {
            Ok(HtsItr{inner: itr, phantom: PhantomData })
        } else {
            Err(hts_err(format!("Error making iterator for file {}", self.name())))
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum HtsRecType {
    Sam, Vcf, Tbx
}

pub struct HtsFile {
    pub(crate) inner: NonNull<htsFile>,
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
    pub fn open<S: AsRef<str>>(name: S, mode: &str) -> io::Result<Self> {
        Self::open_format_(name, mode, None)
    }

    pub fn open_format<S: AsRef<str>>(name: S, mode: &str, fmt: &HtsFormat) -> io::Result<Self> {
        Self::open_format_(name, mode, Some(fmt))
    }

    fn open_format_<S: AsRef<str>>(name: S, mode: &str, fmt: Option<&HtsFormat>) -> io::Result<Self> {
        let name = name.as_ref();
        let mode = get_cstr(mode);
        let fmt = match fmt {
            Some(f) => f.inner.as_ptr(),
            None => null(),
        };
        if let Some(hfile) = NonNull::new(unsafe { hts_open_format(get_cstr(name).as_ptr(), mode.as_ptr(), fmt) }) {
            assert!(!unsafe {hfile.as_ref()}.fn_.is_null());
            let hfile = hfile.cast::<htsFile>();
            Ok(HtsFile {
                inner: hfile,
                phantom: PhantomData,
            })
        } else {
            Err(hts_err(format!("Couldn't open file {}", name)))
        }
    }

    pub fn set_thread_pool(&mut self, tp: &HtsThreadPool) -> io::Result<()> {
        if unsafe{hts_set_thread_pool(self.as_mut(), tp)} == 0 {
            Ok(())
        } else {
            Err(hts_err(format!("Failed to set thread pool for file {}", self.name())))
        }
    }
}

#[derive(Debug)]
pub enum HtsHdr {
    Vcf(VcfHeader),
    Sam(SamHeader),
    Tbx(Tbx),
}

pub trait HtsRead {
    fn read(&mut self, hts: &mut Hts) -> io::Result<bool>;
}

#[repr(C)]
pub struct HtsThreadPool {
    inner: NonNull<hts_tpool>,
    qsize: c_int,
    phantom: PhantomData<hts_tpool>,
}

unsafe impl Send for HtsThreadPool {}

impl Drop for HtsThreadPool {
    fn drop(&mut self) {
        unsafe { hts_tpool_destroy(self.inner.as_ptr() ) }
    }
}

impl HtsThreadPool {
    pub fn new(nthreads: usize) -> Option<Self> {
        NonNull::new(unsafe {hts_tpool_init(nthreads as c_int)})
           .map(|tpool| Self {
                inner: tpool,
                qsize: 0,
                phantom: PhantomData,
            })
    }
}

pub struct HtsIndex {
    inner: NonNull<hts_idx_t>,
    phantom: PhantomData<hts_idx_t>,
}

unsafe impl Send for HtsIndex {}

impl Deref for HtsIndex {
    type Target = hts_idx_t;
    #[inline]
    fn deref(&self) -> &hts_idx_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for HtsIndex {
    #[inline]
    fn deref_mut(&mut self) -> &mut hts_idx_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<hts_idx_t> for HtsIndex {
    #[inline]
    fn as_ref(&self) -> &hts_idx_t { self}
}

impl AsMut<hts_idx_t> for HtsIndex {
    #[inline]
    fn as_mut(&mut self) -> &mut hts_idx_t { self}
}

impl Drop for HtsIndex {
    fn drop(&mut self) {
        unsafe{ hts_idx_destroy(self.as_mut()) }
    }
}

/* impl HtsIndex {

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
                self.idx.idx(),
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
*/

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

pub struct HtsItr{
    inner: NonNull<hts_itr_t>,
    phantom: PhantomData<hts_itr_t>,
}

impl Deref for HtsItr {
    type Target = hts_itr_t;
    #[inline]
    fn deref(&self) -> &hts_itr_t { unsafe{self.inner.as_ref()} }
}

impl DerefMut for HtsItr {
    #[inline]
    fn deref_mut(&mut self) -> &mut hts_itr_t {unsafe{ self.inner.as_mut() }}
}

impl AsRef<hts_itr_t> for HtsItr {
    #[inline]
    fn as_ref(&self) -> &hts_itr_t { self}
}

impl AsMut<hts_itr_t> for HtsItr {
    #[inline]
    fn as_mut(&mut self) -> &mut hts_itr_t { self}
}

impl Drop for HtsItr {
    fn drop(&mut self) {
        unsafe { hts_itr_destroy(self.as_mut()) };
    }
}
