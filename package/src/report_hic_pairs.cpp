#include "diffhic.h"
#include <cctype>

/***********************************************************************
 * Finds the fragment to which each read (or segment thereof) belongs. 
 ***********************************************************************/

class fragment_finder {
public:
	fragment_finder(SEXP, SEXP); // Takes a list of vectors of positions and reference names for those vectors.
	int find_fragment(const int&, const int&, const bool&, const int&) const;
	size_t nchrs() const;
private:
	struct chr_stats {
		chr_stats(const int* s, const int* e, const int& l) : start_ptr(s), end_ptr(e), num(l) {}
		const int* start_ptr;
		const int* end_ptr;
		int num;
	};
	std::deque<chr_stats> pos;
};

fragment_finder::fragment_finder(SEXP starts, SEXP ends) {
	if (!isNewList(starts) || !isNewList(ends)) { throw std::runtime_error("start/end positions should each be a list of integer vectors"); }
	const int nnames=LENGTH(starts);
	if (nnames!=LENGTH(ends)) { throw std::runtime_error("number of names does not correspond to number of start/end position vectors"); }
	
	for (int i=0; i<nnames; ++i) {
		SEXP current1=VECTOR_ELT(starts, i);
		if (!isInteger(current1)) { throw std::runtime_error("start vector should be integer"); }
		SEXP current2=VECTOR_ELT(ends, i);
		if (!isInteger(current2)) { throw std::runtime_error("end vector should be integer"); }
		const int ncuts=LENGTH(current1);
		if (LENGTH(current2)!=ncuts) { throw std::runtime_error("start/end vectors should have the same length"); }
		pos.push_back(chr_stats(INTEGER(current1), INTEGER(current2), ncuts));	
	}
	return;
}

int fragment_finder::find_fragment(const int& c, const int& p, const bool& r, const int& l) const {
	// Binary search to obtain the fragment index with 5' end coordinates. 
	int index=0;
	if (r) {
		int pos5=p+l-1;
		const int* eptr=pos[c].end_ptr;
		index= std::lower_bound(eptr, eptr+pos[c].num, pos5)-eptr;
		if (index==pos[c].num) { 
			warning("read aligned off end of chromosome");
			--index;
		} else if (pos[c].start_ptr[index] > pos5) {
			warning("read aligned into spacer region for a filled-in genome");
			--index;
		}
	} else {
		/* If forward strand, we search relative to the start position of each fragment. If the
 		 * end of the identified fragment is less than the position, the read probably just 
 		 * slipped below the start position of the true fragment and is partially sitting 
 		 * in the spacer for filled-in genomes. We then simply kick up the index.
 		 */
		const int* sptr=pos[c].start_ptr;
		index=std::upper_bound(sptr, sptr+pos[c].num, p)-sptr-1;
		if (pos[c].end_ptr[index] < p) { 
			warning("read starts in spacer region for a filled genome");
			++index; 
		}
	}
	return index;
}

size_t fragment_finder::nchrs() const { return pos.size(); }
	
/***********************************************************************
 * Parses the CIGAR string to extract the alignment length, offset from 5' end of read.
 ***********************************************************************/

void parse_cigar (const char* cigar, int& alen, int& offset, const bool& reverse) {
	alen=0;
	offset=0;
	int numero=0;
    while ((*cigar)!='\0') {
        if (isdigit(*cigar)) {
            numero*=10;
            numero+=int((*cigar)-'0');
        } else {
            switch (*cigar) {
				/* If we've already got some non-clipped values, it means that we're now looking at
			     * the right hand side  of the alignment. Which we don't really need if the read is
			     * forward (as we only want the stuff on the left, as it's equal to the 5' offset).
			     * So we can just return.
				 */
                case 'H': 
					if (alen && !reverse) { return; }					
					offset+=numero;
					break;				
				case 'S':
					if (alen && !reverse) { return; }					
                    break;
                /* We want the length of the alignment relative to the reference, so we ignore insertions
                 * which are only present in the query. We also ignore padded deletions, as they're talking
                 * about something else. We consider deletions from the reference and skipped regions as being  
                 * part of the query. 
                 * 
                 * Note the lack of breaks. It's deliberate, because in all non-clipping cases, we want to check
                 * whether we need to reset the offset (if we're looking for the offset on the right of the 
                 * alignment i.e. the 5' offset on the reverse read). 
                 */
                case 'D': case 'M': case 'X': case '=': case 'N':
                    alen+=numero;
                default:
					if (reverse) { offset=0; }
                    break;
            }
            numero=0;
        }
        ++cigar;
    }
	return;
}

/***********************************************
 * Something to hold segment, pair information.
 ***********************************************/

enum status { ISPET, ISMATE, NEITHER };

struct segment { 
	int offset, alen, fragid, chrid, pos;
	bool reverse;
	bool operator<(const segment& rhs) const { return offset < rhs.offset; }
};

int get_status (const segment& left, const segment& right) {
	if (right.chrid!=left.chrid || right.fragid!=left.fragid || right.reverse==left.reverse) { return NEITHER; }
	const segment& fs=(left.reverse ? right : left);
	const segment& rs=(left.reverse ? left : right);
	if (fs.pos <= rs.pos) {
		if (fs.pos + fs.alen > rs.pos + rs.alen) { return NEITHER; }
		return ISPET; 
	} 
	if (fs.pos < rs.pos+rs.alen) { return NEITHER; }
	return ISMATE;
}

struct valid_pair { 
	int anchor, target, apos, tpos, alen, tlen;
};

/************************
 * Main loop.
 ************************/

SEXP report_hic_pairs (SEXP start_list, SEXP end_list, SEXP pairlen, SEXP chrs, SEXP pos, 
		SEXP flag, SEXP cigar, SEXP mapqual, SEXP chimera_strict, SEXP minqual, SEXP do_dedup) try {
	fragment_finder ff(start_list, end_list);
	const size_t nc=ff.nchrs();

	// Checking input values.
	if (!isInteger(pairlen)) { throw std::runtime_error("length of pairs must be an integer vector"); }
	if (!isInteger(chrs)) { throw std::runtime_error("chromosomes must be an integer vector"); }
	if (!isInteger(pos)) { throw std::runtime_error("positions must be an integer vector"); }
	if (!isInteger(flag)) { throw std::runtime_error("SAM flags must be an integer vector"); }
	if (!isString(cigar)) { throw std::runtime_error("CIGAR strings must be a character vector"); }
	if (!isInteger(mapqual)) { throw std::runtime_error("mapping quality must be an integer vector"); }
	const int nreads=LENGTH(chrs);
	if (LENGTH(pos)!=nreads || LENGTH(flag)!=nreads || LENGTH(cigar)!=nreads || LENGTH(mapqual)!=nreads) { 
		throw std::runtime_error("lengths of vectors of read information are not consistent"); 
	}
	if (!isLogical(chimera_strict) || LENGTH(chimera_strict)!=1) { throw std::runtime_error("chimera removal specification should be a logical scalar"); }
	const int npairs=LENGTH(pairlen);
	if (!isLogical(do_dedup) || LENGTH(do_dedup)!=1) { throw std::runtime_error("duplicate removal specification should be a logical scalar"); }
	if (!isInteger(minqual) || LENGTH(minqual)!=1) { throw std::runtime_error("minimum mapping quality should be an integer scalar"); }

	// Initializing pointers.
	const int* cptr=INTEGER(chrs);
	const int* pptr=INTEGER(pos);
	const int* fptr=INTEGER(flag);
	const int* qptr=INTEGER(mapqual);
	const bool rm_invalid=asLogical(chimera_strict);
	const bool rm_dup=asLogical(do_dedup);
	const int minq=asInteger(minqual);
	const bool rm_min=!ISNA(minq);
	const int * plptr=INTEGER(pairlen);

	// Constructing output containers
	std::deque<std::deque<std::deque<valid_pair> > > collected(nc);
	for (size_t i=0; i<nc; ++i) { collected[i].resize(i+1); }
	std::deque<segment> read1, read2;
	segment current;
	valid_pair curpair;
	int single=0;
	int total=0, dupped=0, filtered=0, mapped=0;
	int dangling=0, selfie=0;
	int total_chim=0, mapped_chim=0, multi_chim=0, inv_chimeras=0;

	// Running through all reads and identifying the interaction they represent.
	int index=0, limit, pindex=0;
	while (index < nreads) {
		read1.clear();
		read2.clear();
		if (pindex==npairs) { throw std::runtime_error("ran out of pairs before running out of reads"); }
		const int& curpl=plptr[pindex];
		++pindex;
		if (curpl==1) { 
			++single;
			++index;
			continue;
		}
		limit=index+curpl;
	    if (limit > nreads) { throw std::runtime_error("ran out of reads before running out of pairs"); }

		// Running through and collecting read segments.
		bool isdup=false, isunmap=false, ischimera=false, skipdup=false, skipunmap=false;
		while (index < limit) {
			const int& curflag=fptr[index];
			current.reverse=(curflag & 0x10);
			current.chrid=cptr[index];
			current.pos=pptr[index];
			parse_cigar(CHAR(STRING_ELT(cigar, index)), current.alen, current.offset, current.reverse);

			// Checking how we should proceed; whether we should bother adding it or not.
			skipdup=(rm_dup && curflag & 0x400);
			skipunmap=(curflag & 0x4 || (rm_min && qptr[index] < minq));
			if (current.offset==0) {
				if (skipdup) { isdup=true; }
				if (skipunmap) { isunmap=true; }
			} else {
				ischimera=true;
			}

			// Checking which deque to put it in, if we're going to keep it.
			if (! skipdup && ! skipunmap) { 
				current.fragid=ff.find_fragment(current.chrid, current.pos, current.reverse, current.alen);
				std::deque<segment>& current_reads=(curflag & 0x40 ? read1 : read2); 
				if (current.offset==0) { 
					current_reads.push_front(current);
				} else {
					current_reads.push_back(current);
				}
			}
			++index;
		}

		/* Adding to statistics, depending on what is going on. In particular, we
		 * skip the read set if the first of either read has any hard 5' clipping. This
 		 * means that it's not truly 5' terminated (e.g. the actual 5' end was unmapped,
 		 * duplicate removed or whatever). Also skipping for any number of other reasons.
 		 */
		++total;
		if (ischimera) { ++total_chim; }
		if (isdup) { ++dupped; }
		if (isunmap) { ++filtered; }
		if (isunmap || isdup || read1.empty() || read2.empty() || read1.front().offset || read2.front().offset) { continue;  }
		++mapped;

		// Determining the type of construct if they have the same ID.
		switch (get_status(read1.front(), read2.front())) {
			case ISPET:
				++dangling;
				continue;
			case ISMATE:
				++selfie;
				continue;
			default:
				break;
		}

		// Pulling out chimera diagnostics.
		if (ischimera) {
			++mapped_chim;
 		   	++multi_chim;	
			bool invalid=false;
			if (read1.size()==1 && read2.size()==1) { 
				--multi_chim;
			} else if (read1.size() > 2 || read2.size() > 2) { 
				invalid=true; 
			} else {
				if (read1.size()==2) { invalid=(ISPET!=get_status(read2[0], read1[1])); }
				if (read2.size()==2) { invalid=(invalid || get_status(read1[0], read2[1])!=ISPET); }
			}
			if (invalid) { 
				++inv_chimeras; 
				if (rm_invalid) { continue; }
			}
		}
		
		// Choosing the anchor segment, and reporting it.
		bool anchor=false;
		if (read1.front().chrid > read2.front().chrid) {
 		   anchor=true;
	   	} else if (read1.front().chrid==read2.front().chrid) { 
			if (read1.front().fragid > read2.front().fragid) {
				anchor=true; 
			} else if (read1.front().fragid == read2.front().fragid) {
				if (read1.front().pos > read2.front().pos) { 
					anchor=true; 
				}
			}
		}
		const segment& anchor_seg=(anchor ? read1.front() : read2.front());
		const segment& target_seg=(anchor ? read2.front() : read1.front());
		
		curpair.anchor=anchor_seg.fragid;
		curpair.target=target_seg.fragid;
		curpair.apos=anchor_seg.pos;
		curpair.alen=anchor_seg.alen; 
		if (anchor_seg.reverse) { curpair.alen*=-1; } 
		curpair.tpos=target_seg.pos;
		curpair.tlen=target_seg.alen;
		if (target_seg.reverse) { curpair.tlen*=-1; }

		if (curpair.alen==0 || curpair.tlen==0) { throw std::runtime_error("alignment lengths of zero should not be present"); } 
		collected[anchor_seg.chrid][target_seg.chrid].push_back(curpair);
	}

	// Checking if all pairs were used up.
	if (pindex!=npairs) { throw std::runtime_error("ran out of reads before running out of pairs"); }

	SEXP total_output=PROTECT(allocVector(VECSXP, 6));
	try { 
		// Checking how many are not (doubly) empty.
		std::deque<std::pair<int, int> > good;
		for (size_t i=0; i<nc; ++i) {
			for (size_t j=0; j<=i; ++j) {
				const std::deque<valid_pair>& curpairs=collected[i][j];
				if (!curpairs.empty()) { good.push_back(std::make_pair(i, j)); }
			}
		}	

		SET_VECTOR_ELT(total_output, 0, allocMatrix(INTSXP, good.size(), 2));
		int* aptr=INTEGER(VECTOR_ELT(total_output, 0));
		int* tptr=aptr+good.size();
		SET_VECTOR_ELT(total_output, 1, allocVector(VECSXP, good.size()));
		SEXP output=VECTOR_ELT(total_output, 1);

		for (size_t i=0; i<good.size(); ++i) {
			aptr[i]=good[i].first+1;
			tptr[i]=good[i].second+1;

			// Filling up those non-empty pairs of chromosomes.
			std::deque<valid_pair>& curpairs=collected[good[i].first][good[i].second];
			SET_VECTOR_ELT(output, i, allocMatrix(INTSXP, curpairs.size(), 6)); 
			int* axptr=INTEGER(VECTOR_ELT(output, i));
			int* txptr=axptr+curpairs.size();
			int* apxptr=txptr+curpairs.size();
			int* tpxptr=apxptr+curpairs.size();
			int* afxptr=tpxptr+curpairs.size();
			int* tfxptr=afxptr+curpairs.size();
			for (size_t k=0; k<curpairs.size(); ++k) {
				axptr[k]=curpairs[k].anchor+1;
				txptr[k]=curpairs[k].target+1;
				apxptr[k]=curpairs[k].apos;
				tpxptr[k]=curpairs[k].tpos;
				afxptr[k]=curpairs[k].alen;
				tfxptr[k]=curpairs[k].tlen;
			}

			// Emptying out the container once we've processed it, to keep memory usage down.
			std::deque<valid_pair>().swap(curpairs);
		}

		// Dumping mapping diagnostics.
		SET_VECTOR_ELT(total_output, 2, allocVector(INTSXP, 4));
		int* dptr=INTEGER(VECTOR_ELT(total_output, 2));
		dptr[0]=total;
		dptr[1]=dupped;
		dptr[2]=filtered;
		dptr[3]=mapped;
	
		// Dumping the number of dangling ends, self-circles.	
		SET_VECTOR_ELT(total_output, 3, allocVector(INTSXP, 2));
		int * siptr=INTEGER(VECTOR_ELT(total_output, 3));
		siptr[0]=dangling;
		siptr[1]=selfie;

		// Dumping the number designated 'single', as there's no pairs.
		SET_VECTOR_ELT(total_output, 4, ScalarInteger(single));

		// Dumping chimeric diagnostics.
		SET_VECTOR_ELT(total_output, 5, allocVector(INTSXP, 4));
		int* cptr=INTEGER(VECTOR_ELT(total_output, 5));
		cptr[0]=total_chim;
		cptr[1]=mapped_chim;
		cptr[2]=multi_chim;
		cptr[3]=inv_chimeras;
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}
	UNPROTECT(1);
	return total_output;
} catch (std::exception& e) {	
	return mkString(e.what());
}

/********************
 * Testing functions.
 *******************/

SEXP test_parse_cigar (SEXP incoming, SEXP reverse) try {
	if (!isString(incoming) || LENGTH(incoming)!=1) { throw std::runtime_error("need one cigar string"); }
	if (!isLogical(reverse) || LENGTH(reverse)!=1) { throw std::runtime_error("need a reverse specifier"); }
	SEXP output=PROTECT(allocVector(INTSXP, 2));
	int* optr=INTEGER(output);
	int& alen=*optr;
	int& offset=*(optr+1);
	parse_cigar(CHAR(STRING_ELT(incoming, 0)), alen, offset, asLogical(reverse));
	UNPROTECT(1);
	return(output);
} catch (std::exception& e) {
	return mkString(e.what());
}

SEXP test_fragment_assign(SEXP starts, SEXP ends, SEXP chrs, SEXP pos, SEXP rev, SEXP len) try {
	fragment_finder ff(starts, ends);
	if (!isInteger(chrs) || !isInteger(pos) || !isLogical(rev) || !isInteger(len)) { throw std::runtime_error("data types are wrong"); }
	const int n=LENGTH(chrs);
	if (n!=LENGTH(pos) || n!=LENGTH(rev) || n!=LENGTH(len)) { throw std::runtime_error("length of data vectors are not consistent"); }
	
	const int* cptr=INTEGER(chrs);
	const int* pptr=INTEGER(pos);
	const int* rptr=LOGICAL(rev);
	const int* lptr=INTEGER(len);

	SEXP output=PROTECT(allocVector(INTSXP, n));
	int *optr=INTEGER(output);
	for (int i=0; i<n; ++i) {
		optr[i]=ff.find_fragment(cptr[i], pptr[i], rptr[i], lptr[i])+1;
	}
	
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

