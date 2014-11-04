/* NINJA OTU Picker: NINJA Is Not Just Another OTU Picker -- filter program
   Knights Lab (www.knightslab.org/ninja)
   This program generates the databases necessary for OTU mapping with bowtie2.
   
   Compilation information (GCC):
   Ascribes to std=gnu and doesn't require C99 support or UNIX intrinsics.
   Flags: -m64 -Ofast -flto ninja_filter.c -fprofile-generate [-fprofile-use]
*/
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#ifndef min 
#define min(a, b) ((a)<=(b) ? (a) : (b)) 
#endif
#define ch(i) *(**(a+i) + depth) 
#define med3(ia, ib, ic) med3func(a, ia, ib, ic, depth)
#define CUTOFF 10
#define MEDCUT 50

/// Utility functions
inline void swap(char ***a, int i, int j) 
	{ char **t = *(a+i); *(a+i) = *(a+j); *(a+j) = t; }

inline int xcmp(str1, str2) register const char *str1, *str2; {
	while (*str1 == *str2++) if (!*str1++) return 0; 
	return (*(const unsigned char *)str1 - *(const unsigned char *)(str2 - 1));
}
inline void vecswap(char ***a, int i, int j, int n) 
	{ while (n-- > 0) swap(a, i++, j++); }

inline int med3func(char ***a, int ia, int ib, int ic, int depth) {
	int va, vb, vc;
	if ((va=ch(ia)) == (vb=ch(ib))) return ia;
	if ((vc=ch(ic)) == va || vc == vb) return ic;
	return va < vb ?
		(vb < vc ? ib : (va < vc ? ic : ia ) ) : 
		(vb > vc ? ib : (va < vc ? ia : ic ) ); 
} 
inline void inssort(char ***a, int n, int depth) {
	int i, j;
	for (i = 1; i < n; i++) for (j = i; j > 0; j--) {
		if (xcmp(**(a+j-1)+depth, **(a+j)+depth) <= 0) break;
		swap(a, j, j-1);
	} 
}  
// 3-way Radix Quicksort (Dobbs-optimized)
void dobbsSrt(char ***a, unsigned n, int depth) {
	if (n < CUTOFF) { inssort(a, n, depth); return; }
	unsigned pl = 0, pm = n >> 1, d;
	int le, lt, gt, ge, r, v, pn = n-1;
	// if large enough, get median of median
	if (n > MEDCUT) {
		d = n >> 3;
		pl = med3(pl, pl+d, pl + (d << 1));
		pm = med3(pm-d, pm, pm+d);
		pn = med3(pn - (d << 1), pn-d, pn);
	}
	pm = med3(pl, pm, pn);
	swap(a, 0, pm);
	v = ch(0); // grab first letter
	for (le = 1; le < n && ch(le) == v; le++);  
	if (le == n) {
		if (v != 0) dobbsSrt(a, n, depth+1);
		return;
	}
	lt = le; gt = ge = n-1;
	// core QS module; partition the data recursively
	for (;;) {
		for ( ; lt <= gt && ch(lt) <= v; lt++)
			if (ch(lt) == v) swap(a, le++, lt);
		for ( ; lt <= gt && ch(gt) >= v; gt--)
			if (ch(gt) == v) swap(a, gt, ge--);
		if (lt > gt) break;
		swap(a, lt++, gt--);
	}
		r = min(le, lt-le); 
		vecswap(a, 0, lt-r, r); 
		r = min(ge-gt, n-ge-1);
		vecswap(a, lt, n-r, r);
		dobbsSrt(a, lt-le, depth);
		if (v != 0) dobbsSrt(a + lt-le, le + n-ge-1, depth+1);
		dobbsSrt(a + n-(ge-gt), ge-gt, depth); 
} 


int cmp(const void *v1, const void *v2) {
	const char *i1 = **(const char ***)v1;
	const char *i2 = **(const char ***)v2;
	return xcmp(i1,i2);
}

int main( int argc, char *argv[] )
{
	if ( argc != 4 ) /* argc should be 4 for correct execution */
    {
        printf( "\nNINJA Is Not Just Another OTU Picker: filter program. Usage:\n");
		printf( "ninja_filter in_reads.fa out_filtered.fa out_sampDB.db\n" );
		printf("\nINPUT PARAMETERS:\n");
		printf( "in_reads.fa: the reads you wish to process\n");
		printf( "\n" "OUTPUT PARAMETERS:\n");
		printf( "out_filtered.fa: the filtered fasta to feed to the aligner (BT2)\n");
		printf( "out_sampDB: bookkeeping DB required by ninja_parse\n");
		return 1;
    }
    //FILE *mfp = fopen( argv[1], "rb" );
	
	char *inmode = "rb", *outmode = "wb"; 
	clock_t start; double cpu_time_used; start = clock();
	// Take in a file with all the user's short reads for later sorting/filtering
	FILE *seqs = fopen(argv[1], inmode);
	if (!seqs) { printf("Error: cannot open input reads fasta.\n"); return 1; }
	FILE *filtered = fopen(argv[2], outmode);
	if (!filtered) {printf("could not open filtered reads output file.\n"); return 1; }
	FILE *sampDBfile = fopen(argv[3], outmode);
	if (!sampDBfile) {printf("could not write to sample DB.\n"); return 1; }
	
	// Find file length
	fseek(seqs, 0, SEEK_END); size_t fsize = ftell(seqs); fseek(seqs, 0, SEEK_SET); 

	char *string = malloc(fsize + 1), *stringP = string - 1; 
	if (!string) { printf("Insufficient memory for short read read-in.\n"); return 1; }
	fread(string, fsize, 1, seqs); //read into string: elems of fsize bytes, 1 elem, using ifp pointer
	fclose(seqs); // close the file
		
	unsigned numEntries = 0, curLen = 1000;  // floating tally
	char *seqBuf = malloc(100000), *seqBufP = seqBuf; // assumption: no read longer than 100KB
	char **seqArr = malloc(curLen * sizeof(char *)), **seqArrP = seqArr,
	*smpBuf = malloc(1024), *smpBufP = smpBuf, // assumption: no sample name > 1KB
	**smpArr = malloc(curLen * sizeof(char *)), **smpArrP = smpArr,
	*thisSeq, *thisSmp;
	if (!smpBuf || !smpArr || !seqArr) 
		{ printf("Cannot allocate post-run memory.\n"); return 1; }
	// loop thru the file again and put things in their place
	unsigned smpCharNum = 0, bufCharNum = 0, inSeq = 0, seqChars = 0;
	stringP = string - 1; // reset read head
	while (*++stringP) { 
		if (*stringP == '>' || *stringP == '-') { continue; } //{}
		else if (*stringP == '\n') { 
			if (!inSeq) { // Header done, prep for sequence
				if (++numEntries == curLen) { // resize all floating arrays
					unsigned offset = seqArrP - seqArr; //, offOtu = DotuP - Dotus;
					//printf("reallocating %u to %u, offset: %u\n", curLen, 2*curLen, offset);
					seqArr = realloc(seqArr, (curLen *= 2) *sizeof(char *));
					smpArr = realloc(smpArr, curLen *sizeof(char *));
					if (!seqArr || !smpArr) 
						{ printf("out of resizing memory.\n"); return 1; }
					seqArrP = seqArr + offset;
					smpArrP = smpArr + offset;
				}
				*smpArrP = malloc(smpCharNum + 1);
#ifdef DEBUG
				if (!*smpArrP) { printf("sample memory depleted!\n"); return 1; }
#endif
				thisSmp = *smpArrP++;
				smpBufP = smpBuf;
				do *thisSmp++ = *smpBufP++; while (--smpCharNum);
				memset(thisSmp, '\0', 1);
				// prepare for inSeq activity
				seqBufP = seqBuf;
				bufCharNum = 0;
			}
			else {
				//memset(seqBufP,'\0',1);
				*seqArrP = malloc(bufCharNum + 1);
#ifdef DEBUG
				if (!*seqArrP) { printf("sequence memory depleted!\n"); return 1; }
#endif
				thisSeq = *seqArrP++;
				seqBufP = seqBuf;
				do *thisSeq++ = *seqBufP++; while (--bufCharNum);
				memset(thisSeq, '\0', 1);
				// prepare for !inSeq activity
				smpBufP = smpBuf;
				smpCharNum = 0;
			}
			inSeq ^= 1;  // swap write target
		}
		else if (inSeq) { 
			*seqBufP++ = *stringP;
			++bufCharNum;
#ifdef DEBUG
			if (bufCharNum > 10000) {printf("sequence buffer overflow\n"); return 1; }
#endif
		}
		else { // in sample header
			if (*stringP == '_') { // fast-forward
				while (*++stringP != '\n'); 
				--stringP; 
				continue;
			} //FF
			*smpBufP++ = *stringP;
			++smpCharNum;
#ifdef DEBUG
			if (smpCharNum > 1024) {printf("sample name buffer overflow\n"); return 1; }
#endif
		}
	}
	seqArr = realloc(seqArr, numEntries * sizeof(char *));
	smpArr = realloc(smpArr, numEntries * sizeof(char *));
	
	
#ifdef PROFILE
	printf("->Short read parse: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	printf("Total short reads: %d\n", numEntries);
	
	// Now sort the sequences and run dedup
	char ***parray = malloc(numEntries * sizeof(char **)), ***pp = parray,
	***smparray = malloc(numEntries * sizeof(char **)), ***smpp = smparray,
	**SmpDD = malloc(numEntries * sizeof(char *)), **ddp = SmpDD;
	
	if (!parray || !smparray || !SmpDD) 
		{ printf("Out of post-memory: parray.\n"); return 1; }
	seqArrP = seqArr; smpArrP = smpArr;
	unsigned nE = numEntries; do { *pp++ = seqArrP++; *smpp++ = smpArrP++; } while (--nE); 
	//qsort(parray, numEntries, sizeof (char *), cmp);
	//qsort(smparray, numEntries, sizeof (char *), cmp);
	dobbsSrt(parray, numEntries, 0);
	dobbsSrt(smparray, numEntries, 0);
	
	// De-dup the sorted sample list and make a parallel counts array
	unsigned copies = 1, dupes = 0; 
	nE = numEntries; smpp = smparray;
	while (--nE) {
		if (!xcmp(**smpp, **(smpp + 1))) ++dupes;
		else *ddp++ = **smpp;
		++smpp;
	}
	unsigned numUniq = numEntries - dupes, nU = numUniq;
	if (ddp - SmpDD <= numUniq) *ddp = **smpp; // endcap
	
	printf("Unique samples: %d\n", numUniq);
	fprintf(sampDBfile, "%u\n", numUniq);
	SmpDD = realloc(SmpDD, numUniq * sizeof(char *)); ddp = SmpDD; 
	nU = numUniq; do fprintf(sampDBfile, "%s\n", *ddp++); while (--nU);
	
	unsigned *Counts = calloc(numUniq, sizeof (unsigned int)), *cpp = Counts;
	if (!Counts) {printf("unable to allocate counts\n"); return 1;}
	
#ifdef PROFILE
	printf("->Short read sample prep: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif

	inline int binSearchS(char *thisS, unsigned range) {
		unsigned binL = 0, binH = range - 1;  
		int binM, cp;
		while (binL <= binH) {
			binM = binL + ((binH - binL) >> 1);
			cp = xcmp(thisS, SmpDD[binM]);
			if (cp > 0) binL = binM + 1;
			else if (!cp) return binM;
			else binH = binM - 1;
		}
		return -1;
	}

	copies = 1; dupes = 0; 
	nE = numEntries; 
	pp = parray; 
	smpBuf = malloc(100000); smpBufP = smpBuf; // assumption: no sample string > 100KB
	unsigned rix = 0;
	int six;
	while (nE--) {
		++*(Counts + (six=binSearchS(*(smpArr + (*pp - seqArr)), numUniq))); // specific count
		if (nE && !xcmp(**pp, **(pp + 1))) { // Dupe! 
			++copies; ++dupes;  // general counts
		} 
		else { // write the last sequence and its number of copies
			smpBufP = smpBuf; 
			six = 0;
			nU = numUniq; cpp = Counts; do {
				if (*cpp) {
					smpBufP += sprintf(smpBufP, "%u:%u:", six, *cpp);
					*cpp = 0;
				}
				++six; ++cpp;
			} while (--nU);
			fprintf(sampDBfile, "%s\n",smpBuf);
			fprintf(filtered, ">%u\n%s\n", rix++, **pp);
			copies = 1;
		}
		++pp;
	}
	printf("Dupes: %d\n", dupes);
	printf("Optimized fasta written.\n");
	// Running is complete
#ifdef PROFILE
	printf("->Read prep and Fasta write: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); 
	start = clock();
#endif
	return 0;
}
