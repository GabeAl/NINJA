#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#define LINELEN 65536
#define SEQPACKS 2048
#define PACKSIZE 32
#define WTYPE uint64_t
#define RSHFT (PACKSIZE*2)-2


#include "qsort.h"
/* typedef struct {
	WTYPE ** addr;
	uint16_t i;
} SortBlock; */
/* inline void SB_qsort(SortBlock *arr, unsigned n) {
	#define SB_LT(a,b) (*(*(a->addr)+(a->i)) < *(*(b->addr)+(b->i)))
	QSORT(SortBlock, arr, n, SB_LT);
} */
typedef struct {
	WTYPE ** addr;
	WTYPE word;
	uint16_t length;
} SortBlock2;

inline void SB2_qsort(SortBlock2 *arr, unsigned n) {
	#define SB2_LT(a,b) (a->word < b->word)
	QSORT(SortBlock2, arr, n, SB2_LT);
}
inline void SB2_qsort2(SortBlock2 *arr, unsigned n) {
	#define SB2_LT2(a,b) ((a->word < b->word) || (a->word == b->word && a->length < b->length))
	QSORT(SortBlock2, arr, n, SB2_LT2);
}

// delete 
//#include "fgets2.c"

WTYPE *C2Xb;
char *X2C = "ACGTNNNNNNNNNNNNNNNN";
char *X2C_RC = "TGCANNNNNNNNNNNNNNNN";
/* 
inline void num2word_order2(WTYPE num, char * word) {
	//char *word = malloc(33);
	//word[32] = 0;
	//unsigned wShf = (chars >> 1) - 2;
	int go = PACKSIZE-1; for (; go > -1; go--) {
		WTYPE temp = num >> RSHFT;
		//printf("%I64u ",temp);
		//if (temp > 3) {puts("ERROR IN WORD."); return;}
		word[go] = X2C[temp];
		num <<= 2;
	}
	//return word;
} */

void num2word(WTYPE num, char * word) {
	//char *word = malloc(33);
	//word[32] = 0;
	//unsigned wShf = (chars >> 1) - 2;
	int go = 0; for (; go < PACKSIZE; go++) {
		WTYPE temp = num >> RSHFT;
		//printf("%I64u ",temp);
		//if (temp > 3) {puts("ERROR IN WORD."); return;}
		word[go] = X2C[temp];
		num <<= 2;
	}
	//return word;
}

void num2wordRC(WTYPE num, char * word) {
	//char *word = malloc(33);
	//word[32] = 0;
	//unsigned wShf = (chars >> 1) - 2;
	int go = PACKSIZE-1; for (; go > -1; go--) {
		WTYPE temp = num >> RSHFT;
		//printf("%I64u ",temp);
		//if (temp > 3) {puts("ERROR IN WORD."); return;}
		word[go] = X2C_RC[temp];
		num <<= 2;
	}
	//return word;
}

char * decodeString(WTYPE * Seq, uint16_t length) {
	//*String = 
	char *newString = malloc(length+1);
	if (!newString) {puts("Error in decodeString."); return 0;}
	newString[length] = 0;
	unsigned clumps = length/PACKSIZE;
	if (PACKSIZE*clumps < length) ++clumps;
	char *word = malloc(33);
	word[32] = 0;
	unsigned z = 0; for (; z < clumps-1; z++) {
		num2word(Seq[z],word);
		memcpy(newString + z*PACKSIZE,word,PACKSIZE);
	}
	num2word(Seq[clumps-1],word);
	memcpy(newString+z*PACKSIZE,word,length % PACKSIZE ?: PACKSIZE); 
	free(word);
	return newString;
}

char * decodeStringX(WTYPE * Seq, uint16_t length, char *word, char *newString) {
	//*String = 
	//char *newString = malloc(length+1);
	//if (!newString) {puts("Error in decodeString."); return 0;}
	newString[length] = 0;
	unsigned clumps = length/PACKSIZE;
	if (PACKSIZE*clumps < length) ++clumps;
	//char *word = malloc(33);
	//word[32] = 0;
	unsigned z = 0; for (; z < clumps-1; z++) {
		num2word(Seq[z],word);
		memcpy(newString + z*PACKSIZE,word,PACKSIZE);
	}
	num2word(Seq[clumps-1],word);
	memcpy(newString+z*PACKSIZE,word,length % PACKSIZE ?: PACKSIZE); 
	//free(word);
	return newString;
}

char * decodeStringRC(WTYPE * Seq, uint16_t length) {
	//*String = 
	char *newString = malloc(length+1);
	if (!newString) {puts("Error in decodeString."); return 0;}
	newString[length] = 0;
	unsigned clumps = length/PACKSIZE;
	if (PACKSIZE*clumps < length) ++clumps;
	char *word = malloc(33);
	word[32] = 0;
	int z = clumps-2; for (; z > -1; z--) {
		num2wordRC(Seq[z],word);
		memcpy(newString + (length) - (z+1)*PACKSIZE,word,PACKSIZE);
	}
	
	num2wordRC(Seq[clumps-1],word);
	//printf("\n%s\n%s\n",word,word+(PACKSIZE - length % PACKSIZE));
	memcpy(newString,word+PACKSIZE-(length % PACKSIZE),(length % PACKSIZE)); 
	free(word);
	return newString;
}

int wcomp(a, b) register const void *a, *b; {
	//return **(WTYPE **)a - **(WTYPE **)b; 
	return **(WTYPE **)a < **(WTYPE **)b ? -1 : **(WTYPE **)a > **(WTYPE **)b; 
}

int comp_wrd(const void *a, const void *b)  {
    if(**((WTYPE **)a) > **((WTYPE **)b))  return(+1);
    if(**((WTYPE **)a) < **((WTYPE **)b))  return(-1);
    return(0);
} 
struct IX2 {
	uint32_t ix1;
	uint32_t ix2;
};
typedef struct {
	WTYPE *W1;
	WTYPE *W2;
} TwoBase;
int baseCmp(Bases, ix) register const void *Bases, *ix; {
	// ix is: size_t
	// Bases is: TwoBase
	return (((TwoBase *)Bases)->W1)[(size_t)ix] - (((TwoBase *)Bases)->W2)[(size_t)ix];
}

/* // slow and easy
int SBCmp(blk1, blk2) register const void *blk1, *blk2; {
	WTYPE **addr1 = ((SortBlock *)blk1)->addr,
		**addr2 = ((SortBlock *)blk2)->addr;
	uint16_t ix1 = ((SortBlock *)blk1)->i,
		ix2 = ((SortBlock *)blk2)->i;
	
	return *(*addr1 + ix1) < *(*addr2 + ix2) ? -1 : *(*addr1 + ix1) > *(*addr2 + ix2);
} */
/* // fast and complicated
int SBCmp2(a,b) register const void *a, *b; {
	return ((*(*((SortBlock*)a)->addr)+(((SortBlock*)a)->i))< 
		*(*(((SortBlock*)b)->addr)+(((SortBlock*)b)->i))) ? -1 : 
		(*(*(((SortBlock*)a)->addr)+(((SortBlock*)a)->i))> 
		*(*(((SortBlock*)b)->addr)+(((SortBlock*)b)->i)));
} */
int SB2Cmp(blk1, blk2) register const void *blk1, *blk2; {
	//WTYPE **addr1 = ((SortBlock *)blk1)->addr,
	//	**addr2 = ((SortBlock *)blk2)->addr;
	//WTYPE w1 = ((SortBlock *)blk1)->word,
	//	w2 = ((SortBlock *)blk2)->word;
	
	//return *(*addr1 + ix1) < *(*addr2 + ix2) ? -1 : *(*addr1 + ix1) > *(*addr2 + ix2);
	return ((SortBlock2 *)blk1)->word < ((SortBlock2 *)blk2)->word ? -1 : 
		((SortBlock2 *)blk1)->word > ((SortBlock2 *)blk2)->word;
}
int SB2Cmp2(blk1, blk2) register const void *blk1, *blk2; {
	return ((SortBlock2 *)blk1)->word < ((SortBlock2 *)blk2)->word ? -1 : 
		((SortBlock2 *)blk1)->word > ((SortBlock2 *)blk2)->word ? 1 :
		((SortBlock2 *)blk1)->length < ((SortBlock2 *)blk2)->length ? -1 :
		((SortBlock2 *)blk1)->length > ((SortBlock2 *)blk2)->length;
}
int SB2Cmp3(blk1, blk2) register const void *blk1, *blk2; {
	if (((SortBlock2 *)blk1)->length != ((SortBlock2 *)blk2)->length) {
	return ((SortBlock2 *)blk1)->word < ((SortBlock2 *)blk2)->word ? -1 : 
		((SortBlock2 *)blk1)->word > ((SortBlock2 *)blk2)->word ? 1 :
		((SortBlock2 *)blk1)->length < ((SortBlock2 *)blk2)->length ? -1 : 1;
		//((SortBlock2 *)blk1)->length > ((SortBlock2 *)blk2)->length;
	}
	return ((SortBlock2 *)blk1)->word < ((SortBlock2 *)blk2)->word ? -1 : 
		((SortBlock2 *)blk1)->word > ((SortBlock2 *)blk2)->word;
}

int SB2Cmp4(blk1, blk2) register const void *blk1, *blk2; {
	//register int w1_lt_w2 = ((SortBlock2 *)blk1)->word < ((SortBlock2 *)blk2)->word;
	//if (w1_lt_w2) return -1;
	if (((SortBlock2 *)blk1)->word < ((SortBlock2 *)blk2)->word) return -1;
	if (((SortBlock2 *)blk1)->word > ((SortBlock2 *)blk2)->word) return 1;
	if (((SortBlock2 *)blk1)->length == ((SortBlock2 *)blk2)->length) return 0;
	if (((SortBlock2 *)blk1)->length < ((SortBlock2 *)blk2)->length) return -1;
	return 1;
}




void superSort(WTYPE ***SeqsP, WTYPE **base, uint16_t *Lengths, 
//inline void superSort(WTYPE **base, uint32_t *SeqsIX, uint16_t *Lengths, 
int depth, size_t beginRange, size_t endRange) {
	// The goal is to read pieces from the sequence pointed to, and swap pointers to sort
	// Step 1: Sort everything in current range by pointer.
	// Copy all refs in this range into a pointer recepticle
	size_t n = endRange - beginRange; // endRange is one after last index
	SortBlock2 *BinPtrs = malloc(n * sizeof(SortBlock2));
	if (!BinPtrs) {puts("Error-MemoryBinPtrs"); return;}
	//depth=1;
	size_t depthSize = (depth+1) * PACKSIZE;
	size_t i = beginRange; for (; i < endRange; ++i) {
		BinPtrs[i-beginRange] = (SortBlock2){SeqsP[i],*(**(SeqsP + i) + depth),
			Lengths[SeqsP[i] - base] <= depthSize ? Lengths[SeqsP[i] - base] : 0};
		
	}
	//SB2_qsort(BinPtrs,n);
	//qsort(BinPtrs, n, sizeof(*BinPtrs), SB2Cmp4);
	SB2_qsort2(BinPtrs,n);
	//FILE *of = fopen("outpiece.txt","wb");
	char *word = calloc(PACKSIZE+1,1);
	if (!word) {puts("error memWORDalloc"); return;}
	
	// Change order of original pointers and free BinPtrs
	for (i=beginRange; i < endRange; ++i) {
		//BinPtrs[i-beginRange] -= depth;
		//num2word(*(*SeqsP[i]+depth),word);
		//fprintf(of,"SeqsP[i]=%I64u (%s) -> ",SeqsP[i],word);
		SeqsP[i] = BinPtrs[i-beginRange].addr; 
		/* if (depth) {
		num2word(*(*SeqsP[i]+depth),word);
		//fprintf(of,"%I64u (%s)\n",SeqsP[i],word);
		printf("\t%I64u: %s\n",i,word);
		} */
	}
	//printf("BASE=%I64u (with first seq addr %I64u).\n",base,*base); 
	//return; 
	/* num2word((*SeqsP[0])[0],word);
	printf("The 0th word in SeqsP now: %I64u=%s\n",(*SeqsP[0])[0],word); */
	
	free(BinPtrs);
	free(word);
	/* FILE *of = fopen("outpiece.txt","wb");
	for (i=beginRange; i < endRange; ++i) {
		fprintf(of,"%I64u\n",*SeqsP[i]);//base[SeqsIX[i]][depth]);
	}
	fclose(of); 
	return; */
	
	// Check for duplicates; for each set, move truncations to top
	WTYPE *curElem = **(SeqsP+beginRange)+depth; // really its address
	size_t lastUniq = beginRange;
	for (i=beginRange + 1; i < endRange; ++i) {
		//if ((*SeqsP[i])[depth] != curElem) {
		if (*(**(SeqsP+i)+depth) != *curElem) {
			//printf("elem %u not equal.\n",i);
			if (i != lastUniq + 1) { //} //printf("Skipping %u\n",i-1);// {  } // // skip range of 1 items
			//else { // loop items here, placing terminal items at top, moving lastU
				size_t z = lastUniq; for (; z < i; ++z) { // float all terminal dupes
					if (Lengths[SeqsP[z] - base] <= depthSize) {
						if (z > lastUniq) {
						//printf("Swapping %I64u with %I64u...\n",z,lastUniq);
						// swap this address for lastUniq++'s address
						WTYPE **temp = SeqsP[z];
						SeqsP[z] = SeqsP[lastUniq];
						SeqsP[lastUniq] = temp;
						}
						//else printf("Sequence %I64u is at its end.\n",z);
						++lastUniq;
					}
				}
				// Spawn a new sort on the remainder
				superSort(SeqsP, base, Lengths, depth+1, lastUniq, i);
				//else ("Error. We have a case where lastUniq >=i.\n");
				
			}
			curElem = **(SeqsP+i)+depth;
			lastUniq = i;
		}
	}
	// endcap
	//--i;
	if (i != lastUniq + 1) { //printf("Skipping ENDCAP %u\n",i-1);
	//else { // loop items here, placing terminal items at top, moving lastU
		size_t z = lastUniq; for (; z < i; ++z) { // float all terminal dupes
			if (Lengths[SeqsP[z] - base] <= depthSize) {
				if (z > lastUniq) {
					//printf("Swapping ENDCAP %I64u with %I64u...\n",z,lastUniq);
					// swap this address for lastUniq++'s address
					WTYPE **temp = SeqsP[z];
					SeqsP[z] = SeqsP[lastUniq];
					SeqsP[lastUniq] = temp;
				}
				//else printf("Sequence ENDCAP %I64u is at its end.\n",z);
				++lastUniq;
			}
		}
		// Spawn a new sort on the remainder
		//if (lastUniq < i) {
			//printf("Spawning new ENDCAP sort (depth %u, range %u to %u)...\n",depth+1,lastUniq,i);
			superSort(SeqsP, base, Lengths, depth+1, lastUniq, i);
		//}
		//else ("Error. We have a case where lastUniq >=i.\n");
		
	}
	
	
}
	
int main () {
	FILE *fp = fopen("seqs.fna", "rb");
	if (fp == NULL) { puts("Invalid input"); return 0; }
	FILE *of = fopen("testout.txt","wb");
	if (of == NULL) { puts("Invalid output"); return 0; }
	
	/* //WSHFT = (((sizeof(WORDTYPE)*4) << 1) - 2); // 62.
	//RSHFT = sizeof(WORDTYPE) * 8 - KMER*2; 
	WTYPE A_ID = (WTYPE)0, C_ID = (WTYPE)1 << RSHFT, 
	G_ID = (WTYPE)2 << RSHFT, T_ID = (WTYPE)3 << RSHFT; 
	C2XbL = calloc(128,sizeof(WTYPE));
	C2XbL['a'] = A_ID; C2XbL['A'] = A_ID; 
	C2XbL['c'] = C_ID; C2XbL['C'] = C_ID; 
	C2XbL['g'] = G_ID; C2XbL['G'] = G_ID; 
	C2XbL['t'] = T_ID; C2XbL['T'] = T_ID; 
	printf("T is %I64u\n",T_ID); */
	char *word = calloc(PACKSIZE+1,1);
	
	C2Xb = calloc(128,sizeof(WTYPE));
	C2Xb['a'] = 0; C2Xb['A'] = 0; 
	C2Xb['c'] = 1; C2Xb['C'] = 1; 
	C2Xb['g'] = 2; C2Xb['G'] = 2;
	C2Xb['t'] = 3; C2Xb['T'] = 3;
	//ctx_t *ctx = init_fgets_sse2 (LINELEN*32); // delete
	//next_t *ne; // delete
	
	/*
	// parse file in 10mb chunks
	char source[MAXBUFLEN + 1];
	source[MAXBUFLEN] = '\0'; // termiNull
	for (;;) {
		size_t newLen = fread(source, 1, MAXBUFLEN, fp);
		if (!newLen) { fputs("File read error", stderr); return 0; }
		//source[++newLen] = '\0'; // termiNull
		
		
		if (n < bufsize) { break; }
	} */
	size_t numElem = 1000, ns=0;
	size_t trim = UINT16_MAX;
	//trim = 30;
	char **Samples = malloc(numElem*sizeof(char *));
	WTYPE **ReadsX = malloc(numElem*sizeof(WTYPE *));
	uint16_t *Sizes = calloc(numElem,sizeof(uint16_t));
	char *line = malloc(LINELEN + 1); // read up to 100k 
	
	while (line = fgets(line,LINELEN,fp)) { // restore
	//while (ne = fgets_sse2(ctx, fp)) { // delete
		//line = ctx->buf; // delete
		if (ns == numElem) {
			//printf("resize beginning\n");
			//while (1);
			numElem *= 2;
			Samples = realloc(Samples,numElem * sizeof(char *));
			ReadsX = realloc(ReadsX, numElem * sizeof(WTYPE *));
			Sizes = realloc(Sizes, numElem*sizeof(uint16_t));
			if (!Samples || !ReadsX || !Sizes) {puts("Error in resize"); return 0;}
			memset(Sizes+numElem/2 + 1,0,(numElem/2-1)*sizeof(uint16_t));
			//size_t z = numElem / 2 + 1; for (; z < numElem; ++z)
			//	Sizes[z] = 0;
		}
		// copy in the sample name up to _ or null minus 1
		char *src = line + 1;
		
		while (*src != '_' && *src != ' ' && *src != '\n') ++src; 
		//printf("\nsrc - line is: %u\n",src-line);

		//size_t count = src - line;
		//Samples[ns] = malloc(LINELEN+1);
		Samples[ns] = malloc(src - line);
		if (!Samples[ns]) {puts("Not enough Samples[ns] mem"); return 0;}
		
		char *dest = Samples[ns]; //char *src = line + 1;
		char *beginSample = line + 1; while (beginSample < src) 
			*dest++ = *beginSample++;
		*dest = 0;
		//memset(dest,0,1);

		// copy in the encoded sequence
		if (!(line = fgets(line,LINELEN,fp))) // restore
		//ne = fgets_sse2(ctx, fp); // delete
		//line = ctx->buf; // delete
		//if (!line) // delete
			{ puts("Error reading file."); return 0; }
		src = line;
		
		//fprintf(of,"Sequence=");
		register size_t length = strlen(src); //, trulen = length;
		if (src[length-1] == '\n') --length; // lop off newline(s)
		if (src[length-1] == '\r') --length; 
		if (trim < length) length = trim;
		
		size_t numPacks = length/PACKSIZE;
		if (numPacks * PACKSIZE < length) ++numPacks;
		
		Sizes[ns] = length; 
		ReadsX[ns] = malloc(numPacks*sizeof(WTYPE));
		if (!ReadsX[ns]) {puts("Bad ReadsX[ns] mem"); return 1; }
		
		WTYPE *thisPack = ReadsX[ns];
		WTYPE clump; int k;
		//fprintf(of,">%s_%u %ubp\n",Samples[ns],ns,Sizes[ns]);
		while (length--) {
			k = 1; clump = C2Xb[*src++];
			while (k < PACKSIZE && length) {
				clump <<= 2u;
				clump += C2Xb[*src++];
				++k; --length;
			}
			if (k != PACKSIZE) clump <<= (2* (PACKSIZE-k));
			*thisPack++ = clump;

			/* num2word(clump,word);
			//if ((thisPack-ReadsX[ns])*PACKSIZE > Sizes[ns]) word[Sizes[ns] % PACKSIZE] = 0;
			if (k < PACKSIZE) word[k] = 0;
			//fprintf(of,"k=%d[%I64u]=%s\n",k,clump,word);
			fprintf(of,"%s",word); */
		}
		//fprintf(of,"\n");
		++ns;
	}
	fclose(fp);
	// Shrink data structures for more memory
	Samples = realloc(Samples,ns * sizeof(char *));
	ReadsX = realloc(ReadsX, ns * sizeof(WTYPE *));
	Sizes = realloc(Sizes, ns * sizeof(uint16_t));
	free(line);
	//printf("Num of sequences: %u\n",ns);
	if (ns > UINT32_MAX) {puts("Too many sequences (>4 bil)."); return 1;}
	printf("max int size=%u/%u\n",sizeof(unsigned),sizeof(uint64_t));
	num2wordRC(ReadsX[0][0],word);
	printf("First number = %I64u\n",ReadsX[0][0]);
	printf("last call=%s on sample: %s, %u elems\n",word,Samples[0], Sizes[0]);
	char *newString = decodeString(ReadsX[0],Sizes[0]);
	char *revString = decodeStringRC(ReadsX[0],Sizes[0]);
	printf("Whole string=\n%s\n",newString);
	printf("Reversed str=\n%s\n",revString);
	//WTYPE word1 = ReadsX[0][0], word2 = ReadsX[1][0];

	//num2word(word1,word); printf("word 1=%s, ",word); 
	//num2word(word2,word); printf(" word2=%s\n",word);
	
	//printf("word 1 bigger than word 2? %d", word1 > word2);
	
	// Create index structure for sequences read (in 32-bit)
	//uint32_t *SeqIX = malloc(sizeof(uint32_t) * ns);
	//uint32_t k = 0; for (; k < ns; ++k) SeqIX[k] = k;
	
	WTYPE ***SeqsP = malloc(sizeof(WTYPE **) * ns);
	size_t k = 0; for (; k < ns; ++k) SeqsP[k] = &ReadsX[k];
	// ...
	//SortBlock a = (SortBlock){ReadsX, 0}, b = (SortBlock){ReadsX+1, 0};
	
	//printf("a's read: %I64u, b's read: %I64u\na comp b? %d\n",**(a.addr),**(b.addr),SBCmp(a,b));
	//num2word(**(a.addr),word);
	//printf("word a=%s, b=",word);
	//num2word(**(b.addr),word);
	//printf("b=%s\n",word);
	
	WTYPE **base = *SeqsP;
	superSort(SeqsP, base, Sizes, 0,0,ns);
	printf("\nDONE SORTING. Printing results to FASTA...\n");
	//char * word = calloc(PACKSIZE+1,1);
	char * string = malloc(UINT16_MAX);
	for (k=0; k < ns; ++k) fprintf(of,">%u\n%s\n",k,
		decodeStringX(*SeqsP[k],Sizes[SeqsP[k]-base],word,string));
	free(SeqsP); // after done using it! */
	
	return 0;
	//uint64_t base = (uint64_t)&ReadsX;
	//printf("First clump is ReadsX[0]=%I64u, or %I64u\n",ReadsX[5][0],*(ReadsX+SeqIX[5])[0]);
	//printf("%d vs %d\n",baseCmp((TwoBase){ReadsX[1],ReadsX[0]},0),(int)(ReadsX[1][0]-ReadsX[0][0]));
	
	//superSort(ReadsX, SeqIX, Sizes, 0, 0, ns);

	return 0;
}
