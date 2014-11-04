/* NINJA OTU Picker: NINJA Is Not Just Another OTU Picker: standalone parser
   OTU table parser, mapper, generator. Uses B-W LF mapping search to
   quickly align each read to the appended DB of all OTU sequences.
   This program assumes the user or companion script has run bowtie2:
   [bowtie2-build-s -o3 ninjaDB.fasta Ninja97] or use included DB
   bowtie2-align-s --no-head --no-unal -o3 -p4 -f reads.fna -x Ninja97 -S align.txt
   [--mp "1,1"--rdg "1,1" --rfg "1,1" --score-min "L,0,-.03" -k 1 --norc --fast]
   
   Compilation information (GCC):
   Ascribes to std=gnu and doesn't require C99 support or UNIX intrinsics.
   Use:     -m64 -Ofast -fwhole-program parse.c -fprofile-generate [-fprofile-use]
   Profile: -m64 -Ofast -D PROFILE -flto parse.c
   Debug:   -m64 -D DEBUG -D PROFILE -ggdb
*/
   
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define ARRSZ 16

inline int ycmp(st1p, st2p) register const char *st1p, *st2p; { 
	while (*st1p && *st2p && *st1p++ == *st2p++); 
	return *--st1p<*--st2p?-1:(*st1p>*st2p); 
}
inline int xcmp(str1, str2) register const char *str1, *str2; {
	while (*str1 == *str2++) if (!*str1++) return 0; 
	return (*(const unsigned char *)str1 - *(const unsigned char *)(str2 - 1));
}

inline int res_unsigned_arr(unsigned **arr, unsigned oldS, unsigned newS) {
	unsigned *newArray = malloc(newS * sizeof(unsigned int)), 
		*nap = newArray, *oap = *arr;
	if (!newArray) { printf("Out of memory.\n"); return 0; }
	int j=-1; while (++j<=oldS) *nap++ = *oap++;
	unsigned * addr = *arr; free(*arr); 
	return (*arr = newArray) - addr;
}
inline int res_str_arr(char ***arr, unsigned oldS, unsigned newS) {
	char **newArray = malloc(newS * sizeof(char *)),
		**nap = newArray, **oap = *arr;
	if (!newArray) { printf("Out of memory.\n"); return 0; }
	int j=-1; while (++j<=oldS) *nap++ = *oap++;
	char ** addr = *arr; free(*arr); 
	return (*arr = newArray) - addr;
}

int cmpS(const void *v1, const void *v2) {
	const char i1 = **(const char **)v1;
	const char i2 = **(const char **)v2;
	return i1<i2?-1:(i1>i2);
}
int cmp(const void *v1, const void *v2) {
	const char *i1 = **(const char ***)v1;
	const char *i2 = **(const char ***)v2;
	return xcmp(i1,i2);
}
int cpx(const void *v1, const void *v2) {
	const int i1 = *(const int *)v1;
	const int i2 = *(const int *)v2;
	return i1<i2?-1:(i1 > i2);
}

int unittest(void)
{
	char *array[ARRSZ] = {"BAD","AS","CATS","FROM","ZEN","XING","WORLDS","UNDER","INVISIBLE","8-balls","BUT","CHOOSING",
		"XEN","MURINE","NAUTILI","PREFERENTIALLY"};
	char **parray[ARRSZ],i;

	for(i = 0; i < ARRSZ;i++) parray[i] = &array[i];
	qsort(parray,ARRSZ, sizeof *parray,cmp);
	puts("Sorting the values:");
	for(i = 0;i < ARRSZ;i++) printf("%s ",array[i]);
	puts("\n");
	for(i = 0;i < ARRSZ;i++)
	printf("value: %s; position: %d\n",*parray[i],
	parray[i]-array);
}

inline unsigned int uWBS(unsigned *ixList, unsigned key, unsigned range) {
	// wide binary search index list for correct OTU 
	unsigned middle, low = 0, high = range;
	while (low <= high) {
		middle = low + ((high - low) >> 1);
		if (key > ixList[middle]) low = middle + 1;
		else if (key < ixList[middle]) high = middle - 1;
		else break; 
	}
	if (ixList[middle] > key) --middle;
	return middle;
}

// file is closed after running
int parse_unsigned_map(FILE * fp, char delim, unsigned ** Col1, unsigned ** Col2) {
	// treat the file to a full read-in and caching
	fseek(fp, 0, SEEK_END); unsigned ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *ixStr = malloc(ixSize + 1); 
	if (!ixStr) { puts("PUM out of memory.\n"); return 0; }
	fread(ixStr, ixSize, 1, fp); fclose(fp);
	unsigned ilines = 0; char *iix = ixStr - sizeof(char), iptr;
	while (iptr = *++iix) iptr != '\n' ?: ++ilines;
	
	// grab values in 'Col1[delim]Col2[\n]' format
	*Col1 = malloc(ilines * sizeof(unsigned int)); 
	*Col2 = malloc(ilines * sizeof(unsigned int));
	if (!*Col1 || !*Col2) { puts("PUM out of memory.\n"); return 0; }
	unsigned *Col1p = *Col1, *Col2p = *Col2;
		 
	int which = 0; // here 0 is the Col1 field and !0 is the Col2 field
	char *buffer = malloc(32), *bufp = buffer;
	iix = ixStr;
	char c; while (c=*iix++) {
		if (c=='\n' || c==delim) { // commit, switch buffer pointer
			memset(bufp,'\0',1);
			which ? (*Col2p++ = atoi(buffer)) : (*Col1p++ = atoi(buffer));
			which ^= 1; bufp = buffer; 
		}
		else *bufp++ = c;
	}
	free(ixStr);
	return ilines;
}
	
// file is closed after running
int parse_string_map(FILE * fp, char delim, unsigned ** Col1, char *** Col2) {
	// treat the file to a full read-in and caching
	fseek(fp, 0, SEEK_END); unsigned ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *ixStr = malloc(ixSize + 1); // *ixStrp = ixStr; 
	if (!ixStr) { puts("PSM out of memory.\n"); return 0; }
	fread(ixStr, ixSize, 1, fp); fclose(fp);
	unsigned ilines = 0; char *iix = ixStr - sizeof(char), iptr;
	while (iptr = *++iix) iptr != '\n' ?: ++ilines;
	
	// grab values in 'Col1[delim]Col2[\n]' format
	*Col1 = malloc(ilines * sizeof(unsigned int)); 
	*Col2 = malloc(ilines * sizeof(char *));
	if (!*Col1 || !*Col2) { puts("PSM out of memory.\n"); return 0; }
	unsigned *Col1p = *Col1; char **Col2p = *Col2;
		
	int which = 0; // here 0 is the Col1 field and !0 is the Col2 field
	char *buffer = malloc(1024), *bufp = buffer;
	iix = ixStr;
	int count = 0;
	char c; while (c=*iix++) {
		if (c=='\n' || c==delim) { // commit, switch buffer pointer
			memset(bufp,'\0',1);
			if (!which) *Col1p++ = atoi(buffer);
			else {
				int amt; *Col2p = malloc(amt = strlen(buffer)+1);
				if (!*Col2p) return 0; 
				strncpy(*Col2p++, buffer, amt);
			}
			which ^= 1; bufp = buffer; 
		}
		else *bufp++ = c;
	}
	free(ixStr);
	return ilines;
}

int main ( int argc, char *argv[] )
{
	//unittest(); return 0;
	clock_t start;
	double cpu_time_used;
	start = clock();
	
    if ( argc != 4 && argc != 5 ) /* argc should be 4-5 for correct execution */
    {
        printf( "\nNINJA Is Not Just Another OTU Picker: standalone parser. Usage:\n");
		printf( "ninja_parse in_aligns.txt in_map.db [in_taxmap.txt] out_otutable.txt\n" );
		printf("\nINPUT PARAMETERS:\n");
		printf( "in_aligns: the bowtie2 (headerless, match-only) short read alignment\n");
		printf( "in_NINJAdb: the (included) sequence index -> OTU database\n");
		printf( "in_taxmap (optional): the (included) sorted OTU -> taxon mapping file\n");
		printf( "\n" "OUTPUT PARAMETERS:\n");
		printf( "out_otutable: name of the new output otu_table file created upon running\n");
		return 1;
    }
    int doTaxmap = argc == 5 ?: 0;
	// We assume argv[n] are filenames to open, in the order specified.
	FILE *ifp = fopen( argv[1], "rb" ), *ifi = fopen( argv[2], "rb"), 
		 *ofp = fopen( argv[argc==4 ? 3 : 4], "wb"), *itx = 0;
	if (doTaxmap) itx = fopen( argv[3], "rb"); // otu-tax map provided
		
	printf("Opened %s for writing\n", argv[3]);
	if ( ifp == 0 || ifi == 0 || ofp == 0 || (doTaxmap && itx == 0))
	{
		fprintf(stderr, "Could not open input or output files.\n");
		printf("usage: ninja_parse in_Map.txt in_NINJAdb.db [in_taxmap.txt] out_otutable.txt");
		return 1;
	}
	
	unsigned *OtuList, *ixList, 
		ilines = parse_unsigned_map(ifi, ',', &ixList, &OtuList);
	char **OtuMap_taxa; unsigned *OtuMap_otus, blines = doTaxmap ? 
		parse_string_map(itx, '\t', &OtuMap_otus, &OtuMap_taxa) : 0;
	if (!ilines || (doTaxmap && !blines)) { printf("Unparsable.\n"); return 1; }
#ifdef PROFILE
	printf("->Time for first 2 funcs: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	fseek(ifp, 0, SEEK_END);
	long fsize = ftell(ifp);
	fseek(ifp, 0, SEEK_SET); 

	char *string = malloc(fsize + 1); // allocate a block of memory (direct bytes)
	if (string == NULL) {
		fprintf(stderr, "Insufficient memory for caching input file.\n");
		return 1;
	}
	fread(string, fsize, 1, ifp); //read into string: elems of fsize bytes, 1 elem, using ifp pointer
	fclose(ifp); // close the file 
#ifdef PROFILE
	printf("->Time for read-in: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	// Give the file a run-thru
	unsigned lines = 0, colons = 0, totCnt = 0, sampCnt = 0, 
		curLen = 2500, linC, startP, lineOTU, SAM = 0, tabs;
	unsigned *Dreps = malloc(curLen * sizeof (unsigned)), *DrepP = Dreps,
		*Dotus = malloc(curLen * sizeof (char *)), *DotuP = Dotus, num;
	char **Dsamps = malloc(curLen * sizeof (char *)), **DsampP = Dsamps;
	char *atoiB = malloc(32 * sizeof(char)), *atoiBp = atoiB,
		*cix = string, *startS; 
	char *head, *dspp, *uhOh = malloc(1024), *uhOhP = uhOh;
	//printf("Size of int on this system: %lu\n",sizeof(__uint128_t));
	//printf("first char is: %c, then: %c\n",*cix, *(cix + 1));
	// Reads either headerless bowtie sam or fasta (respectively, in if/while evals)
	if ((SAM = *cix == ':') || *(cix+++1)==':') do { 
		if (lines < 3) printf("samuel is %d (new tested: %d)", SAM, (*cix == ':'));
		while (*++cix != '\t' && *cix != '\n') { // work through the sample string
			// Single sample name consideration
			startS = cix;
			while (*++cix != ':'); // skip sample name (assumption: at least one char)
			*DsampP = malloc(cix - startS + 1); // prepare to store sample
			head = startS; dspp = *DsampP; 
			do *dspp++ = *head; while (++head < cix);
			//printf("sample here: %s\n", *DsampP);
			++DsampP;
			memset(dspp, '\0', 1);
			
			// Attendant count number parsing
			atoiBp = atoiB;
			while (*++cix != ':') *atoiBp++ = *cix; // read count number
			memset(atoiBp,'\0',1); // commit count string
			num = atoi(atoiB); // convert count string
			totCnt += num; // update total OTU number with new count
			*DrepP++ = num;
			//*DotuP++ = uWBS(ixList, 
			if (curLen == ++sampCnt) { // resize arrays
				unsigned offset = DsampP - Dsamps, offOtu = DotuP - Dotus;
				//printf("reallocating %u to %u, offset: %u\n", curLen, 2*curLen, offset);
				Dsamps = realloc(Dsamps, (curLen *= 2) *sizeof(char *));
				Dreps = realloc(Dreps, curLen *sizeof(unsigned));
				Dotus = realloc(Dotus, curLen *sizeof(unsigned));
				if (!Dsamps || !Dreps || !Dotus) 
					{ printf("out of D-memory.\n"); return 1; }
				DsampP = Dsamps + offset;
				DrepP = Dreps + offset;
				DotuP = Dotus + offOtu;
			}
		}
		if (SAM) { // map OTU directly
			tabs = 0; do while (*++cix != '\t'); while (++tabs != 2);
			startS = cix; while (*++cix != '\t');
			
			uhOhP = uhOh; while (*cix != '\t') *uhOhP++ = *cix++;
			printf("seq is: %s\n", uhOh);
		} 
		//lineOTU = uWBS(
		//startP = sampCnt - (DotuP - Dotus);
		//do *DotuP++ = lineOTU; while (--startP);
		while (*++cix != '\n'); ++cix;
		++lines;
	} while (*cix++); 	
	else do if (*cix == '\n') ++lines; while (*++cix);
	
	printf("Expanded reads: %d (compacted: %d, samps: %d)\n", totCnt, lines, sampCnt);
	//printf("Number of reads: %lu\n", lines >> 1);
	cix = string - 1;
	char **Sample = malloc(lines * sizeof(char *)), **SampleP = Sample;
	unsigned *Otu = calloc(lines, sizeof(unsigned int)), *OtuP = Otu;
	if (!Otu) { printf("Parser module out of memory.\n"); return 1; }
	unsigned count = 0, i; 

#ifdef PROFILE
	printf("->Time for prepass, line count: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	//todo: remove construction/deconstruction in favor of reuse
	while (*++cix) {
		char * idS = malloc(100000), *idSp = idS; // assume < 100KB
		while (*cix != '_' && *cix != '\n') *idSp++ = *cix++;
		memset(idSp,'\0',1);
		// for each sample in the string, create a new entry
		
		*SampleP = malloc(strlen(idS));
		strcpy(*SampleP++,idS);
		free(idS);
		int tabs = 0; while ( tabs < 2 ) *cix++ != '\t' ?: ++tabs;
		char * startS = malloc(32), *startSp = startS;
		while (*cix != '\t') *startSp++ = *cix++;
		memset(startSp,'\0',1);
#ifdef DEBUG
		// Check that the read actually mapped
		if (!ycmp(startS, "*")) { 
			printf("FATAL: Line %lu doesn't map to NINJAdb!\n", count); 
			return 1; 
		}
#endif
		startSp = startS; // reset pointer to re-use field
		while (*++cix != '\t') *startSp++ = *cix;
		memset(startSp,'\0',1);
		
		/* // wide binary search index list for correct OTU 
		unsigned start = atoi(startS);
		unsigned middle, low = 0, high = ilines;
		while (low <= high) {
			middle = low + ((high - low) >> 1);
			if (start > ixList[middle]) low = middle + 1;
			else if (start < ixList[middle]) high = middle - 1;
			else break; 
		}
		if (ixList[middle] > start) --middle; */
		// middle now holds index to the correct OTU ID
		
		//*OtuP++ = OtuList[middle]; //todo: add as many times as there are samples parsed
		*OtuP++ = *(OtuList + uWBS(ixList, atoi(startS), ilines));
		
		free(startS);
		while (*++cix != '\n');
		++count;
	}
	free(string); 
	SampleP = Sample; 
	
#ifdef PROFILE
	printf("->Time for main parse and lookup: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	// sort the sample names by their pointers
	char ***parray = malloc(lines * sizeof(char **)), ***pP = parray; 
	if (!parray) { printf("Insufficient pointer memory.\n"); return 1; }
	linC = lines; do *pP++ = SampleP++; while (--linC);
	qsort(parray,lines, sizeof *parray,cmp);
	
	// rearrange OTUs to match new order
	char **SampleS = malloc(lines * sizeof(char *)), ***pp = parray, **SampleSP = SampleS;
	unsigned *OtuS = malloc(lines * sizeof(unsigned int)), *OtuSP = OtuS;
	linC = lines; do {
		*OtuSP++ = *(Otu + (*pp - Sample)); 
		*SampleSP++ = **pp++;
	} while (--linC); 
	SampleSP = SampleS; OtuSP = OtuS;
	
	// Now that we have an ordered [SampleID, OTU] pairing, we can order OTUs within a 
	// single sample via dynamically resizing piecemeal sort
	unsigned uniqS = 1, thisN = 1, arsz1 = 300,
		*SampBounds = malloc(arsz1 * sizeof (unsigned int)), *SBp = SampBounds;
	char **UniqS = malloc(arsz1 * sizeof (char *)), **USp = UniqS, 
		*thisS = *SampleSP;
	*USp++ = thisS; *SBp++ = 0;
	for (i = 1; i < lines; i++) {
		if (xcmp(*++SampleSP,thisS)) {
			++uniqS;
			qsort(OtuS + i - thisN, thisN, sizeof *OtuS, cpx);
			thisN = 1;
			thisS = *SampleSP;
			*USp++ = thisS; *SBp++ = i;
			if (uniqS > arsz1 - 2) 
				if (!(SBp += res_unsigned_arr(&SampBounds, uniqS, arsz1*=2)) ||
					!(USp += res_str_arr(&UniqS, uniqS, arsz1))) return 1;
		}
		else ++thisN;
	}
	qsort(OtuS + i - thisN,thisN, sizeof *OtuS,cpx); // sort final piece
	*(SampBounds + uniqS) = lines; // append final bounds cap
	
#ifdef PROFILE
	printf("->Piecemeal sort: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	qsort(Otu,lines, sizeof *Otu,cpx); //Otu list now sorted!
	unsigned uniqOtuN = 0, thisOtu = *Otu, arsz = 6000,
		*UniqOtu = malloc(arsz * sizeof(unsigned int)); 
	*UniqOtu = thisOtu; OtuP = Otu;
	linC = lines; while (--linC) if (*OtuP++ != thisOtu) {
		thisOtu = *(OtuP-1);
		*(UniqOtu+ ++uniqOtuN) = thisOtu;
		if (uniqOtuN > arsz - 2) if (!res_unsigned_arr(&UniqOtu, uniqOtuN,arsz*=2))
			return 1;
	}
	++uniqOtuN;
	printf("Unique samples: %d (%d size), unique OTUs: %d (%d size)\n",uniqS,arsz1,uniqOtuN,arsz);
#ifdef PROFILE
	printf("->Time for unique sample calc: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	// build a matrix in memory first. Then write it out in one go.
	thisN = 1; thisS = *SampleS;
	unsigned *Matrix = calloc(uniqS * uniqOtuN, sizeof(unsigned int)),
		curOtu = *UniqOtu, *curBound = &SampBounds[1];
	if (!Matrix) { printf("Cannot allocate matrix.\n"); return 1; }
	char *curSamp = *SampleS;
	int otuIX = -1, sampIX = 0;  
	for (i = 0; i < lines; i++) {
		if (i < *curBound) {
			curOtu = *(OtuS + i);
			while (*(UniqOtu + otuIX) != curOtu) ++otuIX; //if (otuIX >= uniqOtuN) 
			++*(Matrix + (sampIX * uniqS + otuIX));   // increment the count at this OTU in this sample
		}
		else { 
			++sampIX; ++curBound; ++curSamp; //= SampleS[i]; // current sample becomes the next one
			otuIX = 0;
		}
	}
#ifdef PROFILE
	printf("->Time for matrix generation: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	// Write headers in OTU table format
	fprintf(ofp, "#OTU ID");
	for (i=0; i<uniqS; i++) fprintf(ofp,"\t%s",UniqS[i]);
	if (doTaxmap) fprintf(ofp,"\tTaxonomy");
	int step = 0;
	unsigned *lockP = UniqOtu, *stepP = OtuMap_otus - 1,
		*uoP = UniqOtu - 1;
	for (i = 0; i < uniqOtuN; i++) {
		fprintf(ofp, "\n%lu", *++uoP);
		int j; for (j = 0; j < uniqS; j++) fprintf(ofp,"\t%lu", *(Matrix+(j*uniqS + i)));
		if (doTaxmap) {
			// Lock-step algorithm is linear time O(N), not WBS' O(N*logN)
#ifdef DEBUG
			int lstStep = step;
			while (*uoP != OtuMap_otus[step] && step < blines) ++step;
			if (step >= blines) { 
				printf("FATAL: mismatched OTU detected (%lu vs %lu), step: %d\n", UniqOtu[i],*uoP,step);
				fprintf(ofp,"\tNA");
				step = lstStep;
				return 1;
			}
			else fprintf(ofp,"\t%s", *(OtuMap_taxa+step));
#else
			while (*lockP != *++stepP); ++lockP;
			fprintf(ofp,"\t%s", *(OtuMap_taxa + (stepP - OtuMap_otus)));
#endif
		}
	}
#ifdef PROFILE
	printf("->Time for table write: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	printf("Run complete.\n");
	return 0;
}
