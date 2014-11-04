/* NINJA OTU Picker: NINJA Is Not Just Another OTU Picker: parser program
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

int parse_strings(FILE *fp, char *** Strings) {
	fseek(fp, 0, SEEK_END); unsigned ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *Dump = malloc(ixSize + 1); 
	if (!Dump) { puts("PS out of memory.\n"); return 0; }
	fread(Dump, ixSize, 1, fp); fclose(fp);
	unsigned ilines = 0, curLen = 1000; 
	*Strings = malloc(curLen * sizeof(char **));
	char *iix = Dump - 1, *bufP, *Buffer, // = malloc(1000000), *bufP = Buffer,
		**StrP = *Strings;
	while (*++iix) {
		bufP = iix;
		while (*++iix != '\n'); //skips empty lines
		*StrP = malloc(iix - bufP + 1);
		Buffer = *StrP++; 
		do *Buffer++ = *bufP; while (++bufP < iix);
		memset(Buffer,'\0',1);
		if (++ilines == curLen) {
			unsigned offset = StrP - *Strings;
			*Strings = realloc(*Strings, (curLen *= 2) * sizeof(char **));
			if (!*Strings) { printf("out of D-memory.\n"); return 1; }
			StrP = *Strings + offset;
		}
	}
	*Strings = realloc(*Strings, ilines * sizeof( char **));
	return ilines;
}
// file is closed after running
int parse_string_map(FILE * fp, char delim, unsigned ** Col1, char *** Col2) {
	// treat the file to a full read-in and caching
	fseek(fp, 0, SEEK_END); unsigned ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *ixStr = malloc(ixSize + 1); // *ixStrp = ixStr; 
	if (!ixStr) { puts("PSM out of memory.\n"); return 0; }
	fread(ixStr, ixSize, 1, fp); fclose(fp);
	unsigned ilines = 0; char *iix = ixStr - 1, iptr;
	while (iptr = *++iix) iptr != '\n' ?: ++ilines;
	
	// grab values in 'Col1[delim]Col2[\n]' format
	*Col1 = malloc(ilines * sizeof(unsigned int)); 
	*Col2 = malloc(ilines * sizeof(char *));
	if (!*Col1 || !*Col2) { puts("PSM out of memory.\n"); return 0; }
	unsigned *Col1p = *Col1; char **Col2p = *Col2;
		
	int which = 0; // here 0 is the Col1 field and !0 is the Col2 field
	char *buffer = malloc(1024), *bufp = buffer;
	iix = ixStr;
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
	clock_t start;
	start = clock();
	
    if ( argc != 5 && argc != 6 ) /* argc should be 4-5 for correct execution */
    {
        printf( "\nNINJA Is Not Just Another OTU Picker: parser program. Usage:\n");
		printf( "ninja_parse in_sampDB.db in_aligns.txt in_map.db [in_taxmap.txt] out_otutable.txt\n" );
		printf("\nINPUT PARAMETERS:\n");
		printf( "in_sampDB: the sample DB file produced by ninja_filter\n"); 
		printf( "in_aligns: the bowtie2 (headerless, match-only) short read alignment\n");
		printf( "in_NINJAdb: the (included) sequence index -> OTU database\n");
		printf( "in_taxmap (optional): the (included) sorted OTU -> taxon mapping file\n");
		printf( "\n" "OUTPUT PARAMETERS:\n");
		printf( "out_otutable: name of the new output otu_table file created upon running\n");
		return 1;
    }
    int doTaxmap = argc == 6 ?: 0;
	// We assume argv[n] are filenames to open, in the order specified.
	FILE *mfp = fopen( argv[1], "rb" );
	FILE *ifp = fopen( argv[2], "rb" ), *ifi = fopen( argv[3], "rb"), 
		 *ofp = fopen( argv[doTaxmap ? 5 : 4], "wb"), *itx = 0;
	if (doTaxmap) itx = fopen( argv[4], "rb"); // otu-tax map provided
		
	printf("Opened %s for OTU Table writing\n", doTaxmap ? argv[5] : argv[4]);
	if ( !mfp || !ifp || !ifi || !ofp || (doTaxmap && itx == 0))
	{
		fprintf(stderr, "Could not open one or more files.\n");
		printf("Usage: \n");
		printf( "ninja_parse in_sampDB.db in_aligns.txt in_map.db [in_taxmap.txt] out_otutable.txt\n" );
		return 1;
	}
	unsigned *OtuList, *ixList,
		ilines = parse_unsigned_map(ifi, ',', &ixList, &OtuList);
	char **OtuMap_taxa, **SampDBdump; 
	unsigned *OtuMap_otus, blines = doTaxmap ? 
		parse_string_map(itx, '\t', &OtuMap_otus, &OtuMap_taxa) : 0;
	printf("Total OTUs available: %u\n", ilines);
	unsigned slines = parse_strings(mfp, &SampDBdump);
	
	if (!ilines || (doTaxmap && !blines) || !slines) { printf("Unparsable.\n"); return 1; }
#ifdef PROFILE
	printf("->Time for list parse: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	unsigned numSamps = atoi(*SampDBdump++), numReads = slines - numSamps - 1; // max
	char **Seq2samp = SampDBdump + numSamps;
	printf("Number of unique samples: %u, max reads: %u\n", numSamps, numReads);
	//printf("first seq 2 samp: %s\n", Seq2samp[0]);
	
	/// Process SAM file
	fseek(ifp, 0, SEEK_END); size_t fsize = ftell(ifp); fseek(ifp, 0, SEEK_SET); 
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
	
	unsigned *OtuTable = calloc(numSamps * ilines, sizeof (unsigned)), otuIX;
	unsigned tabs, rix, six, cnt, alignPos;
	char *cix = string - 1, *startS, *curSamp; 
	
	while (*++cix) { // != '\t' && *cix != '\n') { // work through the sample string
		startS = cix;
		// lookahead to the read map region
		tabs = 0; do if (*++cix == '\t') ++tabs; while (tabs < 3);
		alignPos = atoi(cix);
		rix = atoi(startS);
		curSamp = *(Seq2samp + rix); // look up rix
		//fprintf(ofp,"Read %u, Align Pos: %u, Sample String: %s\n", rix, alignPos, curSamp);
		otuIX = uWBS(ixList, alignPos, ilines); //*(OtuList + otuIX) is actual otu
		//fprintf(ofp,"-- This read maps to otu %u, (%s)\n", OtuList[otuIX], OtuMap_taxa[uWBS(OtuMap_otus,OtuList[otuIX],blines)]); //*(OtuList + otu));
		do { //startS new scope
			startS = curSamp;
			while (*++curSamp != ':');
			six = atoi(startS);
			startS = ++curSamp; // + 1;
			while (*++curSamp != ':');
			cnt = atoi(startS);
			//fprintf(ofp,"-- Found %u copies in sample id %u (%s).\n", cnt, six, SampDBdump[six] );
			*(OtuTable + otuIX * numSamps + six) += cnt;
		} while (*++curSamp);
		while (*++cix != '\n');
	}
	//return 1;
	// Write headers in OTU table format
	fprintf(ofp, "#OTU ID");
	char **SampP = SampDBdump - 1;
	unsigned *OtuP = OtuList - 1;
	cnt = numSamps; do fprintf(ofp,"\t%s",*++SampP); while (--cnt);
	if (doTaxmap) fprintf(ofp,"\tTaxonomy");
	unsigned i, *row, *rowP;
		
#ifdef PROFILE
	printf("->Time for matrix generation: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
		
		for (i = 0; i < ilines; i++) {
			//screen if any numbers in this row
			row = OtuTable + i*numSamps; rowP = row;
			cnt = numSamps; do if (*rowP++) break; while (--cnt);
			++OtuP;
			if (cnt) {
				fprintf(ofp, "\n%u", *OtuP);
				rowP = row; cnt = numSamps; 
				do fprintf(ofp, "\t%u", *rowP++); while (--cnt);
				if (doTaxmap) 
					fprintf(ofp,"\t%s", 
						*(OtuMap_taxa + uWBS(OtuMap_otus, *(OtuList + i), blines)));
			}
		}

#ifdef PROFILE
	printf("->Time for filtering and table write: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
/*	
#ifdef DEBUG
		// Check that the read actually mapped
		if (!ycmp(startS, "*")) { 
			printf("FATAL: Line %lu doesn't map to NINJAdb!\n", count); 
			return 1; 
		}
#endif
*/	

	printf("Run complete.\n"); 
	return 0;
}
