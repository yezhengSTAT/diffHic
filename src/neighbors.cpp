#include "neighbors.h"

basic::basic(int w, int t, bool i, int x=0) : level(0), width(w), tlen(t), intra(i), exclude(x) {
	if (width < 0 || exclude < 0) { throw std::runtime_error("width values must be non-negative"); }
	if (exclude >= width) { throw std::runtime_error("exclusion width must be less than flank width"); }
}

void basic::restrain () {
	if (left < 0) { left=0; }
	if (intra) {
		if (right > row) { right=row+1; }
	} else if (right > tlen) { right=tlen; } // For intra's, right will hit diagonal; no need to worry about tlen.
	if (left > right) { left=right; } // Avoid negative areas in calling function.
}

/* The bottomright class identifies all bin pairs in a square with sides 'w'.
 * The bin pair of interest lies at the bottomright corner.
 */

bottomright::bottomright(int w, int t, bool i, int x) : basic(w, t, i, x) { level=-w; }

bool bottomright::bump_level () { 
	if (level >= 0) { return false; }
	++level;
	return true;
}

void bottomright::set(int a, int t) {
	row=a+level;
	left=(level < -exclude ? t : t+exclude+1);
	right=t+width+1; 
	restrain();
}

/* The updown class identifies all bin pairs in a vertical line of length 'w*2+1'.
 * The bin pair of interest lies in the centre.
 */

updown::updown(int w, int t, bool i, int x) : basic(w, t, i, x) { level=-w; }
	
bool updown::bump_level() { 
	if (level >= width) { return false; }
	++level;
	if (level==-exclude) { level=exclude+1; }
	return true; 
}

void updown::set(int a, int t) {
	row=a+level;
	left=t;
	right=t+1;
	restrain();
}

/* The leftright class identifies all bin pairs in a horizontal line of length 'w*2+1'.
 * The bin pair of interest lies in the centre.
 * Here, bumping is done to cycle twice over each level; once to get the left side
 * of the line, and again to get the right side.
 */

leftright::leftright(int w, int t, bool i, int x) : basic(w, t, i, x), bumped(false) {} 
	
bool leftright::bump_level() { 
    if (bumped) { 
        return false; 
    } else {
        bumped=true;
        return true;
    }
}

void leftright::set(int a, int t) { 
	row=a;
    if (bumped) {
        left=t+exclude+1;
        right=t+width+1;
    } else {
        left=t-width;
        right=t-exclude;
    }
	restrain(); 
}	

/* `allaround` is split into two rotationally-symmetric shapes; sort of
 * like the L-blocks in tetris, which - when fitted against each other
 * - form a 3x3 square with a hole in the middle. That hole is the 
 * excluded zone in this context, generalized for WxW squares.
 */

allaround::allaround(int w, int t, bool i, int x) : basic(w, t, i, x), bumped(false) { level=-w; } 

bool allaround::bump_level () { 
    if (bumped) { 
        if (level >= width) { return false; }
        ++level;
        return true;
    } else {
        if (level >= width) { // Resetting.
            bumped=true;
            level=-width;
        } else {
           ++level;
        }
        return true;
    }
}

void allaround::set(int a, int t) {
	row=a+level;
    if (bumped) {
    	left=(level > exclude ? t-exclude : t+exclude+1);
        right=t+width+1;
    } else {
        left=t-width;
        right=(level < -exclude ? t+exclude+1 : t-exclude);
    }
    restrain();		
}


