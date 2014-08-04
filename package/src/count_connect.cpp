#include "diffhic.h"
#include "coord.h"

SEXP count_connect(SEXP all, SEXP start, SEXP end, SEXP region, SEXP filter) try {		
	if (!isInteger(start) || !isInteger(end) || !isInteger(region)) { throw std::runtime_error("fragment/indices for first regionion must be integer vectors"); }
	const int ni=LENGTH(start), nr=LENGTH(region)+1; // +1 for 1-based indexing.
	if (LENGTH(end)!=ni) { throw std::runtime_error("start/end index vectors should be the same length"); }

	// Setting up pointers (-1 for 1-based indexing).
	const int * sptr=INTEGER(start)-1,
			* eptr=INTEGER(end)-1,
			* rptr=INTEGER(region)-1;

	// Checking scalars.
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value must be an integer scalar"); }
	const int filtval=asInteger(filter);

	// Setting up other structures, including pointers. We assume it's sorted on R's side.
   	if (!isNewList(all)) { throw std::runtime_error("data on interacting PETs must be contained within a list"); }
	const int nlibs=LENGTH(all);
    std::deque<const int*> aptrs(nlibs), tptrs(nlibs), cptrs(nlibs);
	std::deque<int> nums(nlibs), indices(nlibs);
    std::priority_queue<coord, std::deque<coord>, std::greater<coord> > next;

	for (int i=0; i<nlibs; ++i) {
        SEXP current=VECTOR_ELT(all, i);
        if (!isNewList(current) || LENGTH(current)!=3) { 
			throw std::runtime_error("interactions must be supplied as a data.frame with anchor.id, target.id and counts"); }

        for (int j=0; j<3; ++j) {
            SEXP current_col=VECTOR_ELT(current, j);
            if (!isInteger(current_col)) { throw std::runtime_error("interaction data must be in integer format"); }
            int* ptr=INTEGER(current_col);
            switch (j) {
                case 0: 
					aptrs[i]=ptr; 
					nums[i]=LENGTH(current_col);
					break;
                case 1: tptrs[i]=ptr; break;
                case 2: cptrs[i]=ptr; break;
                default: break;
            }
		}
        // Populating the priority queue.
		if (nums[i]) { next.push(coord(aptrs[i][0], tptrs[i][0], i)); }
	}
	
	// Running through all libraries.
	std::deque<int> counts;
	typedef std::pair<int, int> combo;
	std::map<combo, std::pair<int, int> > bins;
	std::map<combo, std::pair<int, int> >::iterator itb;
	std::deque<int> curcounts(nlibs);

	int curab, curtb, curlib;
	combo temp;
	int counter=0;
	while (!next.empty()) {
		curab=next.top().anchor;
		curtb=next.top().target;
		++counter;

		// This speeds things up by collecting entries with the same fragment indices.
		do {
			curlib=next.top().library;
			int& libdex=indices[curlib];
			curcounts[curlib]+=cptrs[curlib][libdex];
			next.pop();
            if ((++libdex) < nums[curlib]) {
				next.push(coord(aptrs[curlib][libdex], tptrs[curlib][libdex], curlib));
			} 
		} while (!next.empty() && next.top().anchor==curab && next.top().target==curtb);

		/* Allocating counts to every pair of ranges containing these fragments. This
 		 * shouldn't be too sadistic if there aren't too many overlapping ranges.
 		 */
		if (curab>ni) { throw std::runtime_error("invalid anchor index for supplied fragments"); } // 1-based indexing, so '>' is right.
		const int& s1x=sptr[curab];
		const int& e1x=eptr[curab];
		if (curtb>ni) { throw std::runtime_error("invalid target index for supplied fragments"); }
		const int& s2x=sptr[curtb];
		const int& e2x=eptr[curtb];
	
		if (s1x!=e1x && s2x!=e2x) { 
			if (s1x <= 0 || s2x <= 0 || e1x > nr || e2x > nr) { throw std::runtime_error("invalid start/endpoints for region indices"); }
			for (int x1=s1x; x1<e1x; ++x1) {
				for (int x2=s2x; x2<e2x; ++x2) { 
					if (rptr[x1] > rptr[x2]) {
						temp.first=rptr[x1];
						temp.second=rptr[x2];
					} else {
						// Avoid redundant naming, when looking within the same ranges.
						temp.first=rptr[x2];
						temp.second=rptr[x1];
					}
					
					itb=bins.lower_bound(temp);
					if (itb==bins.end() || bins.key_comp()(temp, itb->first)) {
						itb = bins.insert(itb, std::make_pair(temp, std::make_pair(counts.size(), counter)));
						counts.resize(counts.size()+nlibs);
					} else if ((itb->second).second==counter) { 
						/* The 'counter' avoids adding the same range twice to a particular pair. Regions
						 * can be irregularly sized and spaced, e.g., nested, so it's not possible to set up 
						 * a general rule to avoid redundant counting.
						 */
						continue;
					}
					(itb->second).second=counter;
					const int& index=(itb->second).first;
					for (int lib=0; lib<nlibs; ++lib) { counts[index+lib] += curcounts[lib]; }
				}
			}
		}

//		/* Okay, now flipping and adding the ones where anchor goes to region2, and target goes to region1.
//		 * This is necessary if there are two ranges, in which case we don't know how they'll react.
//		 */
//		if (!issame) {
//			if (curtb>ni1) { throw std::runtime_error("invalid target index for supplied fragments"); } 
//			const int& s1y=s1ptr[curtb];
//			const int& e1y=e1ptr[curtb];
//			if (curab>ni2) { throw std::runtime_error("invalid anchor index for supplied fragments"); }
//			const int& s2y=s2ptr[curab];
//			const int& e2y=e2ptr[curab];
//
//			if (s1y!=e1y && s2x!=e2y) {
//				// Running through, and skipping if there are clashes. 
//				for (int y1=s1y; y1<e1y; ++y1) {
//					temp.first=r1ptr[y1];
//					for (int y2=s2y; y2<e2y; ++y2) { 
//						temp.second=r2ptr[y2];
//						
//						itb=bins.lower_bound(temp);
//						if (itb==bins.end() || bins.key_comp()(temp, itb->first)) {
//							itb = bins.insert(itb, std::make_pair(temp, std::make_pair(counts.size(), counter)));
//							counts.resize(counts.size()+nlibs);
//						} else if ((itb->second).second==counter) { // Don't add the same range twice. 
//							continue;
//						}
//						(itb->second).second=counter;
//						const int& index=(itb->second).first;
//						for (int lib=0; lib<nlibs; ++lib) { counts[index+lib] += curcounts[lib]; }
//					}
//				}
//			}
//		}
		
		// Resetting.
		for (int i=0; i<nlibs; ++i) { curcounts[i]=0; }
	}

	// Assessing how many combinations are above threshold.
	int index=0, countsum=0, ncombos=0;
	std::deque<bool> isokay;
	while (index < counts.size()) { 
		for (int i=0; i<nlibs; ++i) { 
			countsum += counts[index];
			++index; 
		}
		if (countsum>=filtval) { 
			isokay.push_back(true);
			++ncombos;
		} else {
			isokay.push_back(false);
		}
		countsum=0;
	}

	// Returning all count combinations underneath the threshold.
	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, ncombos));
		int* aoptr=INTEGER(VECTOR_ELT(output, 0));
		SET_VECTOR_ELT(output, 1, allocVector(INTSXP, ncombos));
		int* toptr=INTEGER(VECTOR_ELT(output, 1));
		SET_VECTOR_ELT(output, 2, allocMatrix(INTSXP, ncombos, nlibs));
		std::deque<int*> coptrs(nlibs);
		for (int i=0; i<nlibs; ++i) {
			if (i==0) { coptrs[i]=INTEGER(VECTOR_ELT(output, 2)); }
			else { coptrs[i]=coptrs[i-1]+ncombos; }
		}	
		
		// Iterating across and filling both the matrix and the components.
		int odex=0;
		for (itb=bins.begin(); itb!=bins.end(); ++itb) {
			if (isokay[(itb->second).first/nlibs]) { 
				aoptr[odex]=(itb->first).first;
				toptr[odex]=(itb->first).second;
				const int& index=(itb->second).first;
				for (int i=0; i<nlibs; ++i) { coptrs[i][odex]=counts[index+i]; }
				++odex;
			}
		}
	} catch (std::exception& e) { 
		UNPROTECT(1);
		throw;
	}

	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}
