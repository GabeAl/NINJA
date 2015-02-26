/* NINJA EMBALMER: Where Alignment Goes to Die (it's only an acronym in Arabic)
   Extensive Multi-ref Burst ALignment by Matrix-Exhaustive Recursion
   Knights Lab (www.knightslab.org/ninja)
   
   Compilation information (GCC):
   Ascribes to std=gnu and doesn't require C99 support or UNIX intrinsics.
   Flags: -m64 -Ofast -flto ninja_EMBALM.c
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

unsigned parse_fasta(FILE *fp, char **header, char **seq) {
	if (!fp) return 0;
	fseek(fp, 0, SEEK_END); unsigned sz = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *str = malloc(sz + 1), *strP; 
	if (!str) { puts("fasta parser out of memory.\n"); return 0; }
	fread(str, sz, 1, fp); fclose(fp);
	strP = str; do; while (*++strP != '\n'); // '\r' can get in here on Windows
	memset(strP,'\0',1); memset(str+sz,'\0',1);
	*header = str; *seq = ++strP;
	return sz - (strP - str);
}

#define LOGGER() { \
	FILE *log = fopen("logE.txt", "ab" ); \
	if (1) { \
		for (i = 0; i < qs; i++) {  \
			for (j = 0; j < rs; j++) fprintf(log,"%d\t",scoreMatrix[i*rs + j]); \
			fprintf(log,"\n"); \
		} \
		fprintf(log,"\n"); \
		for (i = 0; i < qs; i++) {  \
			for (j = 0; j < rs; j++) fprintf(log,"%d\t",trackMatrix[i*rs + j]); \
			fprintf(log,"\n"); \
		} \
		fprintf(log,"................\n"); \
	} \
}

float alignRQ(char *ref, char *query, unsigned rs, unsigned qs, int ends, int print) {
	const unsigned char DIAG = 1, UP = 2, LEFT = 3;
	//int gap = -1, miss=-1, match=1; // hard-coded for id generation
	++qs; ++rs;
	int *scoreMatrix = calloc(qs * rs, sizeof(int));
	unsigned char *trackMatrix = malloc(qs * rs * sizeof(int));
	int i, j; 
	for (i = 1; i < qs; i++) scoreMatrix[i*rs] = -i; //gap * i; // starting gaps in query penalized
	if (ends) for (i = 1; i < rs; i++) scoreMatrix[i] = -i;
	for (i = 0; i < rs; i++) trackMatrix[i] = LEFT; // first row starts pointing left
	for (i = 0; i < qs; i++) trackMatrix[i*rs] = UP; // first col starts pointing up
	int scoreD, scoreU, scoreL, dir, maxScore;
	for (i = 1; i < qs; i++) for (j = 1; j < rs; j++) {
		scoreD = scoreMatrix[(i-1)*rs + j-1]; // diag [format is [row*rowsize + col]]
		scoreU = scoreMatrix[(i-1)*rs + j];
		scoreL = scoreMatrix[i*rs + j-1];
		if (query[i-1] == ref[j-1] || query[i-1] == 'N' || ref[j-1] == 'N') ++scoreD;// += match; 
		else --scoreD; // += miss;
		--scoreU; --scoreL; //+= gap; // gap penalty with every shift up/left
		maxScore = scoreD; dir = DIAG; // Decide max score and attendant shift
		if (scoreL > maxScore) { maxScore = scoreL; dir = LEFT; }
		if (scoreU > maxScore) { maxScore = scoreU; dir = UP; }
		scoreMatrix[i*rs + j] = maxScore; // best direction's score
		trackMatrix[i*rs + j] = dir; // record best direction
	}
	
	// find the max score relative to the query
	int maxSize = rs * 2 + 1, botIX = (qs-1)*rs, maxI = INT_MIN;
	unsigned T = 0, maxIX=0; // alignment size thus far
	for (i=0; i < rs; i++) if (scoreMatrix[botIX + i] > maxI) 
		{ maxI = scoreMatrix[botIX + i]; maxIX = i; }
	if (print) { // generate a full graphical alignment
		char *seqA = malloc(maxSize), *seqB = malloc(maxSize);
		seqA[0] = ref[ends ? (rs - 1) : maxIX];
		seqB[0] = query[qs-1];
		i = qs - 1; j = ends ? (rs - 1) : maxIX; 
		
		while (i || (ends && j)) { // count down until 0,0... or until query if !ends
			if (trackMatrix[i*rs + j] == DIAG) {
				seqA[T] = ref[j---1];
				seqB[T] = query[i---1];
			}
			else if (trackMatrix[i*rs + j] == UP) {
				seqA[T] = '-';
				seqB[T] = query[i---1];
			}
			else { // trackMatrix == LEFT
				seqA[T] = ref[j---1];
				seqB[T] = '-';
			}
			++T;
		}
		char *alQuery = malloc(T+1), *alRef = malloc(T+1), *bar = malloc(T+1);
		for (i = 0; i < T; i++) {
			alRef[i] = seqA[T-i-1];
			alQuery[i] = seqB[T-i-1];
			bar[i] = alQuery[i] == alRef[i] ? '|' : ' ';
		}
		memset(alQuery+T,'\0',1); memset(alRef+T,'\0',1); memset(bar+T,'\0',1);
		printf("%s\n%s\n%s\n", alRef, bar, alQuery);
	}
	else { // perform a fast traceback
		i = qs - 1; j = ends ? (rs - 1) : maxIX; 
		while (i || (ends && j)) {
			if (trackMatrix[i*rs + j] == DIAG) { --i; --j; }
			else if (trackMatrix[i*rs + j] == UP) --i;
			else --j;
			++T;
		}
	}
#ifdef LOG_MATRICES
	LOGGER();
#endif
	return 1.f - (T - (ends ? scoreMatrix[qs*rs-1] : maxI))*.5/T;
}

int main( int argc, char *argv[] )
{
	FILE *ifp1, *ifp2; // init file pointers (if applicable)
#ifdef LOG_MATRICES
	FILE *log = fopen("logE.txt", "wb" ); fclose(log); // clear existing log
#endif
	char *refHead, *refSeq, *qHead, *qSeq;
	unsigned refLen, qLen;
	
	
	if (argc < 3 || argc > 4) { // Neither use-case satisfied
		puts("Usage: ref.fa query.fa\n"); 
		puts("OR: REFERENCESEQ QUERYSEQ output_appended_scores.txt\n");
		
		// Do a toy example
		refLen = 11; qLen = 7;
		refSeq = "AGGATACGTCA\0";
		qSeq = "AGATCGA\0"; 
		printf("Match id: %f\n\n", alignRQ(refSeq, qSeq, refLen, qLen, 0, 1));
		printf("Align id: %f\n__________________\n", alignRQ(refSeq, qSeq, refLen, qLen, 1, 1));
		char *ref2 = "AGGATCGTCA\0";
		char *q2 = "GGATACGA\0";
		printf("Match id: %f\n\n", alignRQ(q2, ref2, 8, 10, 0, 1));
		printf("Align id: %f\n", alignRQ(q2, ref2, 8, 10, 1, 1));
	}
	if (argc == 3) { // This is the standard Filename1 Filename2 case
		ifp1 = fopen( argv[1], "rb" );
		ifp2 = fopen( argv[2], "rb" );
		refLen = parse_fasta(ifp1, &refHead, &refSeq);
		qLen = parse_fasta(ifp2, &qHead, &qSeq);
		if (!refLen || !qLen) {puts("Invalid input file(s)!\n"); return 1; }
		printf("Match id: %f\n\n", alignRQ(refSeq, qSeq, refLen, qLen, 0, 1));
		printf("Align id: %f\n__________________\n", alignRQ(refSeq, qSeq, refLen, qLen, 1, 1));
	}
	else if (argc == 4) { // Assume REF QUERY outfile.txt use-case
		FILE *outs = fopen(argv[3], "ab");
		if (!outs) { puts("Invalid output file. Is it open?\n"); return 1; }
		float score = alignRQ(argv[1], argv[2], strlen(argv[1]), strlen(argv[2]), 0, 0);
		fprintf(outs, "%f\n", score);
		printf("%f\n", score);
	}
	
	return 0;
}
