use std::{
    io,
    marker::PhantomData,
    ptr::NonNull,
    ops::{Deref, DerefMut},
    convert::{AsRef, AsMut, TryFrom},
    ffi::CStr,
    collections::HashMap,
};

pub mod lib;
pub use lib::*;

use super::{
    from_cstr, c_to_cstr, get_cstr, htsFile, hts_err, kstring_t, Hts, HtsFile, HtsPos, HtsHdr,
    HtsRead, HtsWrite, HtsItr, BGZF, hts_itr_next, HtsFileDesc
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
    pub fn new(mode: &str) -> io::Result<Self> {
        match NonNull::new(unsafe { bcf_hdr_init(get_cstr(mode).as_ptr()) }) {
            None => Err(hts_err("Couldn't create VCF/BCF header".to_string())),
            Some(hdr) => Ok(VcfHeader { inner: hdr , phantom: PhantomData}),
        }
    }

    pub fn dup(&self) -> io::Result<Self> {
        match NonNull::new(unsafe { bcf_hdr_dup(self.as_ref()) }) {
            None => Err(hts_err("Couldn't duplicate VCF/BCF header".to_string())),
            Some(hdr) => Ok(VcfHeader { inner: hdr , phantom: PhantomData}),
        }
    }

    pub fn read(fp: &mut HtsFile) -> io::Result<Self> {
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

    pub fn seq_names(&self) -> Vec<&str> {
        let mut n_seq: c_int = 0;
        let p = unsafe{bcf_hdr_seqnames(self.as_ref(), &mut n_seq as *mut c_int)};
        if p.is_null() {
            Vec::new()
        } else {
            let mut v = Vec::with_capacity(n_seq as usize);
            for i in 0..n_seq {
                let c_str: &CStr = unsafe { CStr::from_ptr(*p.offset(i as isize)) };
                let str_slice: &str = c_str.to_str().unwrap();
                v.push(str_slice);
            }
            unsafe {libc::free(p as *mut c_void)};
            v
        }
    }

    pub fn seq_lengths(&self) -> Vec<usize> {
        let n = self.n_ref() as usize;
        (0..n).map(|i| self.id2len(i)).collect()
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

    pub fn fmt<'a, T: BcfHeaderVar<'a, T>>(&self, tag: &str) -> io::Result<BcfTag<T>> { BcfTag::fmt(self, tag) }

    pub fn info<'a, T: BcfHeaderVar<'a, T>>(&self, tag: &str) -> io::Result<BcfTag<T>> { BcfTag::info(self, tag) }

    pub fn flt(&self, tag: &str) -> io::Result<BcfTag<BcfFlag>> { BcfTag::flt(self, tag) }

    pub fn flt_set(&self, tags: &[&str]) -> io::Result<FilterSet> { FilterSet::new(self, tags) }
}

impl Clone for VcfHeader {
    fn clone(&self) -> Self {
        self.dup().expect("Error duplicating VCF/BCF header")
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
unsafe impl Sync for BcfRec {}

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

impl HtsWrite for BcfRec {
    fn write(&mut self, hts: &mut Hts) -> io::Result<()> {
        let (fp, hdr) = hts.hts_file_and_header();
        let hts_file = fp.as_mut();
        let res = if let Some(HtsHdr::Vcf(hd)) = hdr {
            let i = unsafe { bcf_write(hts_file, hd.as_mut(), self.as_mut()) };
            match i {
                0 => Ok(()),
                _ => Err(hts_err("Error writing VCF/BCF record".to_string())),
            }
        } else {
            Err(hts_err(format!("Appropriate header not found for file {}", fp.name())))
        };
        res
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

    pub fn n_allele(&self) -> u16 { self.as_ref().n_allele() }
    pub fn n_info(&self) -> u16 { self.as_ref().n_info() }
    pub fn n_fmt(&self) -> u8 { self.as_ref().n_fmt() }
    pub fn n_flt(&mut self) -> c_int {
        if (self.unpacked as usize) & BCF_UN_FLT == 0 { self.unpack(BCF_UN_FLT)}
        self.d.n_flt
    }

    pub fn passed(&mut self) -> bool { self.n_flt() > 0 }

    pub fn qual(&self) -> f32 { self.as_ref().qual }
    pub fn n_sample(&self) -> u32 { self.as_ref().n_sample() }

    pub fn format(&mut self, hdr: &HtsHdr, s: &mut kstring_t) -> io::Result<()> {
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

    pub fn get_alleles<'a>(&'a mut self) -> Vec<&'a str> {
        self.unpack(BCF_UN_STR);
        let n_all = self.n_allele() as usize;
        let mut v = Vec::with_capacity(n_all);
        let all = &self.d.alleles;
        if all.is_null() { panic!("BCF allele desc is null")}
        for i in 0..n_all {	v.push(from_cstr::<'a>(unsafe{*all.add(i)}))}
        v
    }
    pub fn get_alleles_by<'a, F>(&'a mut self, mut f: F)
    where
        F: FnMut(&'a [u8])
    {
        self.unpack(BCF_UN_STR);
        let n_all = self.n_allele() as usize;
        let all = &self.as_ref().d.alleles;
        if all.is_null() { panic!("BCF allele desc is null")}

        for i in 0..n_all {
            f(c_to_cstr( unsafe { *all.add(i) } ).to_bytes())
         }
    }

    pub fn get_fmt<'a, T: BcfHeaderVar<'a, T>>(&'a mut self, htag: &BcfTag<T>) -> Option<BcfVecIter<'a, T>> { htag.get_fmt(self) }

    pub fn get_info<'a, T: BcfHeaderVar<'a, T>>(&'a mut self, htag: &BcfTag<T>) -> Option<BcfValIter<'a, T>> { htag.get_info(self) }

    pub fn get_flt(&mut self, htag: &BcfTag<BcfFlag>) -> bool { htag.get_flt(self) }

    pub fn flt_set_any(&mut self, fset: &FilterSet) -> bool { fset.any(self) }

    pub fn flt_set_all(&mut self, fset: &FilterSet) -> bool { fset.all(self) }

}

impl Drop for BcfRec {
    fn drop(&mut self) {
        unsafe { bcf_destroy(self.inner.as_ptr()) }
    }
}

pub struct FilterSet {
    hset: HashMap<c_int, u64>,
    mask: u64,
}

impl FilterSet {
    pub fn new(hdr: &VcfHeader, tags: &[&str]) -> io::Result<Self> {
        if tags.len() > 64 {
            return Err(hts_err("Too many tags: Maximum tags in a FilterSet is 64".to_string()))
        }
        let tag_ids: Vec<_> = tags.iter().enumerate().filter_map(|(ix, s)| {
             hdr.id2int(BCF_DT_ID as usize, s)
               .and_then(|id| {
                   let info = unsafe {(*(*hdr.as_ref().id[BCF_DT_ID as usize].add(id)).val).info[BcfHeaderLine::Flt as usize] };
                   if (info & 0xf) == 0xf { None } else { Some((ix, id)) }
               })
        }).collect();
        if tag_ids.len() < tags.len() {
            let mut missing_tags = Vec::with_capacity(tags.len() - tag_ids.len());
            let mut ix = 0;
            for (i, _) in tag_ids.iter() {
                for s in &tags[ix..*i] { missing_tags.push(*s)}
                ix = *i + 1;
            }
            for s in &tags[ix..] { missing_tags.push(*s)}
            Err(hts_err(format!("Not all filter tags were found in header.  Missing tags: {:?}", missing_tags)))
        } else {
            let hset: HashMap<c_int, u64> = tag_ids.iter().copied().enumerate().map(|(ix, (_, i))| (i as c_int, 1u64 << ix)).collect();
            let mask = (1u64 << tags.len()) - 1;
            Ok(Self{hset, mask})
        }
    }

    // Check if any of the filters in FilterSet are present
    pub fn any(&self, rec: &mut BcfRec) -> bool {
        if (rec.unpacked as usize) & BCF_UN_FLT == 0 { rec.unpack(BCF_UN_FLT)}
        assert!(rec.d.n_flt >= 0);
        let n_flt = rec.d.n_flt as usize;
        if n_flt == 0 && self.hset.contains_key(&0) { true }
        else {
            let p = rec.d.flt;
            assert!(!p.is_null());
            unsafe { std::slice::from_raw_parts(p, n_flt) }
               .iter().any(|i| self.hset.contains_key(i))
        }
    }

    // Check if all of the filters in FilterSet are present
    pub fn all(&self, rec: &mut BcfRec) -> bool {
        if (rec.unpacked as usize) & BCF_UN_FLT == 0 { rec.unpack(BCF_UN_FLT)}
        assert!(rec.d.n_flt >= 0);
        let n_flt = rec.d.n_flt as usize;
        if n_flt == 0 { self.mask == 1 }
        else {
            let p = rec.d.flt;
            assert!(!p.is_null());
            self.mask == unsafe { std::slice::from_raw_parts(p, n_flt) }
               .iter().fold(0, |m, i| if let Some(x) = self.hset.get(i) { m | x } else { m })
        }
    }
}

pub struct BcfTag<T> {
    tag_hl: BcfHeaderLine,
    tag_id: isize,
    n_samples: usize,
    phantom: PhantomData<T>,
}

impl <'a, T: BcfHeaderVar<'a, T>> BcfTag<T> {
    pub fn new(hdr: &VcfHeader, mut tag: &str, tag_hl: BcfHeaderLine) -> io::Result<BcfTag<T>> {
        if tag_hl > BcfHeaderLine::Fmt {
            return Err(hts_err("Bad Header line type used for BcfHeaderVar".to_string()))
        }
        if tag_hl == BcfHeaderLine::Flt && tag == "." { tag = "PASS" }
        let hlt = tag_hl as usize;
        let tag_id = hdr.id2int(BCF_DT_ID as usize, tag)
           .ok_or_else(|| hts_err(format!("Unknown header tag {}", tag)))? as isize;
        let info = unsafe {(*(*hdr.as_ref().id[BCF_DT_ID as usize].offset(tag_id)).val).info[hlt] };
        if (info & 0xf) == 0xf {
            return Err(hts_err(format!("Unknown header tag {} of type {:?}", tag, tag_hl)))
        }
        let ty = if tag_hl == BcfHeaderLine::Flt { BcfHeaderType::Flag } else {
            BcfHeaderType::try_from((info >> 4) & 0xf).expect("Illegal header var type")
        };
        if tag == "GT" {
            if ty != BcfHeaderType::Str { return Err(hts_err(format!("Incorrect header variable type for GT tag {}", tag))) }
            if T::hdr_type() != BcfHeaderType::Int { return Err(hts_err(format!("Incorrect variable type for GT tag {}", tag))) }
        } else if ty != T::hdr_type() {
            return Err(hts_err(format!("Incorrect variable type for header tag {}", tag)))
        }
        let n_samples = hdr.n_samples();
        assert!(n_samples >= 0);

        Ok(Self{tag_hl, tag_id, n_samples: n_samples as usize, phantom: PhantomData})
    }

    pub fn fmt(hdr: &VcfHeader, tag: &str) -> io::Result<BcfTag<T>> { Self::new(hdr, tag, BcfHeaderLine::Fmt) }

    pub fn info(hdr: &VcfHeader, tag: &str) -> io::Result<BcfTag<T>> { Self::new(hdr, tag, BcfHeaderLine::Info) }

    pub fn flt(hdr: &VcfHeader, tag: &str) -> io::Result<BcfTag<T>> { Self::new(hdr, tag, BcfHeaderLine::Flt) }

    pub fn get_info(&self, rec: &'a mut BcfRec) -> Option<BcfValIter<'a, T>> {
        assert_eq!(self.tag_hl, BcfHeaderLine::Info, "Wrong tag type - expected Info tag");
        if (rec.unpacked as usize) & BCF_UN_INFO == 0 { rec.unpack(BCF_UN_INFO)}
        let n_info = rec.n_info() as usize;
        if n_info > 0 {
            let p = rec.d.info;
            assert!(!p.is_null());
            unsafe { std::slice::from_raw_parts(p, n_info) }
               .iter().find(|i| i.key == (self.tag_id as c_int))
               .and_then(|i| {
                   NonNull::new(i.vptr).map(|p| {
                       assert!(i.len > 0 && i.type_size() > 0, "BCF vector length or size is zero!");
                       let (n, size) = if i.vtype == BCF_BT_CHAR {
                           (1, i.len as usize)
                       } else {
                           (i.len as usize, i.type_size())
                       };
                       BcfValIter{inner: BcfGenIter{p, n, size, phantom: PhantomData}, missing: false}
                   })})

        } else { None }
    }

    pub fn get_fmt(&self, rec: &'a mut BcfRec) -> Option<BcfVecIter<'a, T>> {
        assert_eq!(self.tag_hl, BcfHeaderLine::Fmt, "Wrong tag type - expected Format tag");
        if (rec.unpacked as usize) & BCF_UN_FMT == 0 { rec.unpack(BCF_UN_FMT)}
        let n_fmt = rec.n_fmt() as usize;
        if n_fmt > 0 {
            let p = rec.d.fmt;
            assert!(!p.is_null());
            unsafe { std::slice::from_raw_parts(p, n_fmt) }
               .iter().find(|f| f.id == (self.tag_id as c_int))
               .and_then(|f| {
                   NonNull::new(f.p).map(|p| {
                       assert!(f.n > 0 && f.size > 0 && f.type_size() > 0, "BCF vector length or size is zero!");
                       let (n, type_size) = if f.vtype == BCF_BT_CHAR {
                           (1, f.n as usize)
                       } else {
                           (f.n as usize, f.type_size())
                       };
                       BcfVecIter{ inner: BcfGenIter{p, n, size: type_size, phantom: PhantomData}, size: f.size as usize, n_vec: self.n_samples }
                   })})

        } else { None }
    }

    pub fn get_flt(&self, rec: &'a mut BcfRec) -> bool {
        assert_eq!(self.tag_hl, BcfHeaderLine::Flt, "Wrong tag type - expected Filter tag");
        if (rec.unpacked as usize) & BCF_UN_FLT == 0 { rec.unpack(BCF_UN_FLT)}
        assert!(rec.d.n_flt >= 0);
        let n_flt = rec.d.n_flt as usize;
        if self.tag_id == 0 && n_flt == 0 { true }
        else {
            let p = rec.d.flt;
            assert!(!p.is_null());
            unsafe { std::slice::from_raw_parts(p, n_flt) }
               .iter().any(|i|  *i == (self.tag_id as c_int))
        }
    }
}

pub struct BcfGenIter<'a, T> {
    p: NonNull<u8>,
    n: usize,
    size: usize,
    phantom: PhantomData<&'a T>,
}

pub struct BcfVecIter<'a, T> {
    inner: BcfGenIter<'a, T>,
    size: usize,
    n_vec: usize,
}

impl <'a, T>Iterator for BcfVecIter<'a, T>
where T: BcfHeaderVar<'a, T>,
{
    type Item = BcfValIter<'a, T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.n_vec == 0 {
            None
        } else {
            let iter = BcfValIter{ inner: BcfGenIter{p: self.inner.p, n: self.inner.n, size: self.inner.size, phantom: PhantomData}, missing: false};
            self.inner.p = unsafe { NonNull::new_unchecked(self.inner.p.as_ptr().add(self.size)) };
            self.n_vec -= 1;
            Some(iter)
        }
    }
}

pub struct BcfValIter<'a, T> {
    inner: BcfGenIter<'a, T>,
    missing: bool,
}

impl <'a, T>Iterator for BcfValIter<'a, T>
    where T: BcfHeaderVar<'a, T>,
{
    type Item = Option<T::Item>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.inner.n == 0 {
            None
        } else {
            self.inner.n -= 1;
            Some(
                if self.missing { None } else {
                    let p = unsafe {std::slice::from_raw_parts(self.inner.p.as_ref(), self.inner.size)};
                    let x = match T::try_parse(p).expect("BCF parse error") {
                        BcfOpt::Some(c) => Some(c),
                        BcfOpt::Missing => None,
                        BcfOpt::EndOfVec => {
                            self.missing = true;
                            None
                        },
                    };
                    self.inner.p = unsafe { NonNull::new_unchecked(self.inner.p.as_ptr().add(self.inner.size)) };
                    x
                }
            )
        }
    }
}

pub fn hts_get_vcf_header(hts: &Hts) -> Option<&VcfHeader> {
    hts.header().and_then(|h| {
        match h {
            HtsHdr::Vcf(vh) => Some(vh),
            HtsHdr::Sam(_) => None,
        }
    })
}