/**
 ** Copyright (c) 2011-2014 Illumina, Inc.
 **
 ** This file is part of the BEETL software package,
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 ** Lightweight BWT Construction for Very Large String Collections.
 ** Proceedings of CPM 2011, pp.219-231
 **
 **/

#include "TransposeFasta.hh"
#include <zlib.h>
#include "kseq.h"

#include "Filename.hh"
#include "SeqReader.hh"
#include "Tools.hh"
#include "libzoo/util/Logger.hh"
#include "libzoo/util/TemporaryFilesManager.hh"

#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <numeric>
#include <functional>
#include <math.h>

//TODO spostare tutto in TransposeFasta.hh oppure parameters
using namespace std;
KSEQ_INIT(gzFile, gzread);
unsigned myTotalLen;
unsigned maxSeqLen;
unsigned minSeqLen;
vector <int> keepOrder;
unsigned nTotalSeq;
unsigned int* keepOrderRLO;


// struct used to establish permutation of sequences by length
struct seqInfo{
public:
    int len;   //sequence length
    unsigned nSeq; //original sequence's position in input file (0 is first)
};
// short function to compare two seqInfo by length, from longest to shortest
bool compareSequence (seqInfo s1, seqInfo s2){
   return s1.len > s2.len;
}

double calcolaVarianza(vector<double> v){
    //varianza e media
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double var = (sq_sum/v.size() - mean*mean);
    double stdev = std::sqrt(var);
    cerr<<"Mean "<<mean<<" Varianza "<<var<<" StdDev "<<stdev<<endl;
    exit(EXIT_SUCCESS);
}

struct charElem{
    char ch;
    unsigned int oldPosition;
};

bool compareCharElem(charElem ch1, charElem ch2){
    return ch1.ch < ch2.ch;
}

//RLO ordering
struct group{
public:
    unsigned int base; //from which position the group starts
    unsigned int bound; //limit position of the group
    vector<charElem> newCharacters;
   group(){
        this->base = 0;
        this->bound = 0;
    }
};

TransposeFasta::TransposeFasta()
    : pReader_( NULL )
    , cycleNum_( 0 )
    , processQualities_( false )
{
    for ( int i( 0 ); i < 256; i++ ) freq[i] = 0;
}

unsigned TransposeFasta::findInfoSeq(const string &input) {

    for (dataTypedimAlpha z = 0 ; z < SIZE_ALPHA-1; z++)
        freq[z]=0;
    freq[SIZE_ALPHA-1]=0;
    freq[(unsigned int)(TERMINATE_CHAR)]=1;

   // vector <double> v;
    vector <seqInfo> infoVector;
    string bufChar;
    int nSeq = 0; //number of sequences found
    long int lenTot = 0; //total characters read
    int sizeAlpha = 0; //different alphabet symbols found
    minSeqLen = UINT_MAX;
    maxSeqLen = 0;

    //calculate file name
    std::string cutNameInput = (input.substr(input.find_last_of("/\\") + 1));
    std::string::size_type const p(cutNameInput.find_last_of('.'));
    std::string file_without_extension = cutNameInput.substr(0, p);

    gzFile fp;
    kseq_t *seq;
    long int l;
    fp = gzopen(input.c_str(), "r");
    if(fp == Z_NULL){
        printf("Error while opening file %s\n", input.c_str());
        exit(EXIT_FAILURE);
    }
    seq = kseq_init(fp);
    int charsBuffered = 0;
    int seqRead = 0;
    long previousLen = -1;

    //while there are new sequences to read
    while ((l = kseq_read(seq)) >= 0) {
        //detect different length among sequences
        if( !differentLenDetected_ && (previousLen >= 0) ){
            if(l != previousLen)
                differentLenDetected_ = true;
        }
        previousLen = l;
        nSeq++;
        for (unsigned int z = 0; z < l; z++)
            freq[(unsigned int) (seq->seq.s[z])] = 1;
        unsigned int actualLen = l;
        lenTot += actualLen;
        //find max
        if (actualLen > maxSeqLen)
            maxSeqLen = actualLen;
        //find min
        if (actualLen < minSeqLen)
            minSeqLen = actualLen;
        seqInfo addSeq;
        addSeq.nSeq = nSeq - 1;
        addSeq.len = actualLen;
        infoVector.push_back(addSeq);
      //  v.push_back(l);
    }
    //calcolaVarianza(v);

    if( ordering == 2 ) {
        sort(infoVector.begin(), infoVector.end(), compareSequence);
        assert( minSeqLen == infoVector[nSeq - 1].len);
        assert( maxSeqLen == infoVector[0].len);
    }

    //sizeAlpha
    for (dataTypedimAlpha i = 0; i < SIZE_ALPHA - 1; ++i)
        if (freq[i] > 0) {
            sizeAlpha++;
        }
    if (freq[SIZE_ALPHA - 1] > 0) {
        sizeAlpha++;
    }

    //Creating info file //TODO metti nella cartella random
    std::stringstream infoFileName;
    //string cutNameInput = input.substr(0, strlen(input.c_str())-4);

    infoFileName << TEMP_DIR << "/" << file_without_extension.c_str() << ".info";
    FILE* outputInfo = fopen( infoFileName.str().c_str(),"w" );
    if(!outputInfo) {
        printf("Error while opening file\n");
        exit(EXIT_FAILURE);
    }
    //Writing into info file
    std::stringstream infoData;
    infoData << "MinLengthRead=" << minSeqLen << "\n" << "MaxlengthRead="
             << maxSeqLen << "\n" << "lenTot=" << lenTot << "\n" << "nSeq="
             << nSeq << "\n" << "sizeAlpha=" << sizeAlpha;
    string data = infoData.str();
    size_t written = fwrite(data.c_str(), sizeof(char), data.length(), outputInfo);
    assert(written == data.length());
    fclose(outputInfo);

    nTotalSeq = nSeq;
if( ordering == 2) {
    //Print permutation by len
    //Creating permutation file
    std::stringstream permFileName;
    permFileName << TEMP_DIR << "/"
                << file_without_extension.c_str()
                << "_len_permutation.txt";
    FILE *outputPerm = fopen(permFileName.str().c_str(), "w");
    if (!outputPerm) {
        printf("Error while opening file\n");
        exit(EXIT_FAILURE);
    }
    for (long int i = 0; i < nTotalSeq; i++)
        fprintf(outputPerm, "%d\n", infoVector[i].nSeq);

    fflush(outputPerm);
    fclose(outputPerm);
}

    myTotalLen = lenTot;
    kseq_destroy(seq);
    gzclose(fp);
    return maxSeqLen;
}

bool TransposeFasta::init( SeqReaderFile *pReader, const string &input, const bool processQualities)
{
    differentLenDetected_ = false;
    //printf("\nFile name %s", input.c_str());
    pReader_ = pReader;
#if (ACCEPT_DIFFERENT_LEN == 1)
    cycleNum_ = findInfoSeq(input); //TODO remove keepOrder and use file instead
    outputFiles_.resize( cycleNum_ );
    buf_.resize( cycleNum_, vector<uchar>( BUFFERSIZE ) );
    if(differentLenDetected_) {
        for (int i = 0; i < cycleNum_; i++)
            for (int j = 0; j < BUFFERSIZE; j++)
                buf_[i][j] = DUMMY_CHAR;
    }
    processQualities_ = processQualities;
#else
    cycleNum_ = pReader->length();
    outputFiles_.resize( pReader->length() );
    buf_.resize( pReader->length(), vector<uchar>( BUFFERSIZE ) );
    processQualities_ = processQualities;
#endif

    cerr << "Constructing TransposeFasta, found max read length of "
         << cycleNum_ << " and min read length of " << minSeqLen << endl;

    if ( processQualities_ && pReader_->thisQual()[0] == '\0' )
    {
        // If the first entry of the file (which can be fastq or any other format (raw/fasta/etc)) doesn't contain any quality info
        // , deactivate qualities processing
        processQualities_ = false;
    }

    return differentLenDetected_;
}

TransposeFasta::~TransposeFasta()
{

}

string TransposeFasta::sortInputFileByLen(const string &input, const string &output, bool generatedFilesAreTemporary) {

    vector < FILE * > lenOutput;
    bool isFirstTime[cycleNum_];
    bool fileExists[cycleNum_];
    for (unsigned i = 0; i < cycleNum_; ++i) {
        isFirstTime[i] = true;
        fileExists[i] = false;
    }
    lenOutput.resize(cycleNum_);

    //Open input file for reading
    gzFile fp;
    kseq_t *seq;
    fp = gzopen(input.c_str(), "r");
    seq = kseq_init(fp);
    string strToWrite;
    int len = 0;

    //Read one sequence at a time and write in the corresponding len file
    while ((len = kseq_read(seq)) >= 0) {
      //  cerr<<len<<endl;
        unsigned index = len - 1;
        //build the string to write
        if (!isFirstTime[index]) {
            strToWrite.append("\n");
        }
        //append sequence id
        strToWrite.append(">");
        if( seq->name.l )
            strToWrite.append(seq->name.s);
        //append comment
        if( seq->comment.l ) {
            strToWrite.append(" ");
            strToWrite.append(seq->comment.s);
        }
        strToWrite.append("\n");
        strToWrite.append(seq->seq.s);
        //append quality
        if (processQualities_ && (seq->qual.l) ) {
            strToWrite.append("\n");
            strToWrite.append(seq->qual.s);
        }

      //check if file exists, otherwise create and set existence
      if(!fileExists[index]) {
          Filename fn(output, index, "");
          lenOutput[index] = fopen(fn, "w");
          if (lenOutput[index] == NULL) {
              cerr << "Error: couldn't open len output file " << fn << endl;
              if (index > 0) {
                  cerr
                          << "  You may have reached the maximum number of "
                             "opened files (see `ulimit -n`) or the maximum number of "
                             "files allowed in one directory"
                          << endl;
                  exit(-1);
              }
          }
          if (generatedFilesAreTemporary)
              TemporaryFilesManager::get().addFilename(fn);
          fileExists[index] = true;
      }
        //write the sequence
        FILE *whereToWrite = lenOutput[index];
        size_t written_bytes = fwrite(strToWrite.c_str(),
                                      sizeof(char),
                                      strToWrite.size(),
                                      whereToWrite);
        assert(written_bytes == strToWrite.size());
        isFirstTime[index] = false;
        strToWrite.clear();
    }

    assert(lenOutput.size() > 1);

    //Create name of new input file
    std::string cutNameInput = (input.substr(input.find_last_of("/\\") + 1));

    //concat files, build cat command
    std::ostringstream command;
    command << "cat ";
    Filename cm(output, "");
    // Order is from longest to shortest
    for(unsigned i = maxSeqLen-1; i > 0; --i) {
        if(fileExists[i]) {
            command << cm << i << " ";
            //print a new line so I can concat with new line
            fprintf(lenOutput[i], "\n");
            fflush(lenOutput[i]);
        }
    } //case 0 outside because is unsigned
    if(fileExists[0]) {
        command << cm << 0 << " ";
        fflush(lenOutput[0]);
    }
    command << "> ";
    command << TEMP_DIR << "/" << cutNameInput.c_str();
   // cerr<<command.str().c_str()<<endl;
    int res = system(command.str().c_str());

    //remove all len files
    for(unsigned i = minSeqLen-1; i < maxSeqLen; ++i){
        if(fileExists[i]) {
            fclose(lenOutput[i]);
            Filename rm(output, i, "");
            remove(rm);
        }
    }

    kseq_destroy(seq);
    gzclose(fp);
    return cutNameInput.c_str();
}

/* Creates cyc files possibly inserting DUMMY_CHAR to fill length difference between strings */
bool TransposeFasta::convert( const string &input, const string &output, bool generatedFilesAreTemporary )
{
    vector<vector<uchar> > bufQual;
    vector<FILE *> outputFilesQual;
    if ( processQualities_ )
    {
        bufQual.resize( cycleNum_, vector<uchar>( BUFFERSIZE ) );
        outputFilesQual.resize( cycleNum_ );
    }

    //TO DO
    lengthRead = cycleNum_;
    //The distribution of characters is useful
    //for alpha[256] -->Corresponding between the alphabet, the piles and tableOcc
    //and to know sizeAlpha
    //We supposed that the symbols in the input file are the following
    freq[int( terminatorChar )] = 1;
    freq[int( 'A' )] = 1;
    freq[int( 'C' )] = 1;
    freq[int( 'G' )] = 1;
    freq[int( 'N' )] = 1;
    freq[int( 'T' )] = 1;
    //GIOVANNA: ADDED THE SYMBOL Z IN THE ALPHABET, SO sizeAlpha = alphabetSize
#ifdef USE_EXTRA_CHARACTER_Z
    freq[int( 'Z' )] = 1;
#endif

    // create output files
    for ( SequenceLength i = 0; i < cycleNum_; i++ )
    {
        Filename fn( output, i, "" );
        outputFiles_[i] = fopen( fn, "w" );
        if ( outputFiles_[i] == NULL )
        {
            cerr << "Error: couldn't open output file " << fn << endl;
            if ( i > 0 )
            {
                cerr << "  You may have reached the maximum number of opened files (see `ulimit -n`) or the maximum number of files allowed in one directory, as we create one file per cycle (and a second one if qualities are present)" << endl;
                exit ( -1 );
            }
        }
        if ( generatedFilesAreTemporary )
            TemporaryFilesManager::get().addFilename( fn );
        if ( processQualities_ )
        {
            Filename fnQual( output + "qual.", i, "" );
            outputFilesQual[i] = fopen( fnQual, "w" );
            if ( outputFilesQual[i] == NULL )
            {
                cerr << "Error: couldn't open output file " << fnQual << endl;
                if ( i > 0 )
                {
                    cerr << "  You may have reached the maximum number of opened files (see `ulimit -n`) or the maximum number of files allowed in one directory, as we create one file per cycle (and a second one if qualities are present)" << endl;
                    exit ( -1 );
                }
            }
            if ( generatedFilesAreTemporary )
                TemporaryFilesManager::get().addFilename( fnQual );
        }
    }

//fill buffer with '#' to manage different length sequences
#if (ACCEPT_DIFFERENT_LEN == 1)
    if(differentLenDetected_) {
        for (int i = 0; i < cycleNum_; i++) {
            for (int j = 0; j < BUFFERSIZE; j++)
                buf_[i][j] = DUMMY_CHAR;
        }
    }
    if ( processQualities_ ){
        for (int i = 0; i < cycleNum_; i++) {
            for (int j = 0; j < BUFFERSIZE; j++)
                bufQual[i][j] = DUMMY_CHAR;
        }
    }
#endif

    // looping through the input file, add the characters to the buffer, print buffer when it's full
    //    unsigned int num_read = 0;
    unsigned int num_write = 0;
    unsigned int charsBuffered = 0;

    lengthTexts = 0;
    nSeq = 0;

    int len = 0;
    gzFile fp;
    kseq_t *seq;
    fp = gzopen(input.c_str(), "r");
    seq = kseq_init(fp);

    while ( (len = kseq_read(seq)) >= 0 ){
        //cerr << "current line : " << buf << endl;
        lengthTexts+=len;
        if ( charsBuffered == BUFFERSIZE )
        {
            // write buffers to the files, clear buffers
#pragma omp parallel for num_threads(4)
            for ( SequenceLength i = 0; i < cycleNum_; i++ )
            {
                //cerr << "writing to " << i << " : " << buf_[i] << endl;
                size_t num_write_bases = fwrite ( buf_[i].data(), sizeof( char ), charsBuffered, outputFiles_[i] );
                checkIfEqual( num_write_bases, charsBuffered ); // we should always read/write the same number of characters
                if ( processQualities_ )
                {
                    size_t num_write_qual = fwrite ( bufQual[i].data(), sizeof( char ), charsBuffered, outputFilesQual[i] );
                    checkIfEqual( num_write_bases, num_write_qual );
                }
            }
            //lengthTexts += ( num_write * cycleNum_ );
// reset buffer with #
#if (ACCEPT_DIFFERENT_LEN == 1)
            if(differentLenDetected_) {
                for (int i = 0; i < cycleNum_; i++) {
                    for (int j = 0; j < BUFFERSIZE; j++)
                        buf_[i][j] = DUMMY_CHAR;
                }
            }
            if ( processQualities_ ){
                for (int i = 0; i < cycleNum_; i++) {
                    for (int j = 0; j < BUFFERSIZE; j++)
                        bufQual[i][j] = DUMMY_CHAR;
                }
            }
#endif

            charsBuffered = 0;
        }

#if (ALIGN == 0) //Right
        //Align right side to compute preprocessing RLO
            int index = cycleNum_-1;
            for ( int i = len-1; i >= 0; --i ){
                buf_[index][charsBuffered] = seq->seq.s[i];

                if (processQualities_ ){
                    bufQual[index][charsBuffered] = seq->qual.s[i];
                }
                index--;
            }
            if(index >= 0)
                buf_[index][charsBuffered] = TERMINATE_CHAR;
#else
        //else align left side
        SequenceLength i;
        for ( i = 0; i < len; ++i ) {
            buf_[i][charsBuffered] = seq->seq.s[i];

            if (processQualities_) {
                bufQual[i][charsBuffered] = seq->qual.s[i];
            }
        }
        if(i <= cycleNum_-1)
            buf_[i][charsBuffered] = TERMINATE_CHAR;

#endif

        // increase the counter of chars buffered
        charsBuffered++;
        nSeq++;

    }

    // write the rest
#pragma omp parallel for num_threads(4) //aggiunto #pragma
    for ( SequenceLength i = 0; i < cycleNum_; i++ )
    {
        num_write = fwrite ( buf_[i].data(), sizeof( uchar ), charsBuffered, outputFiles_[i] );
        lengthTexts += num_write;
        if ( processQualities_ )
        {
            size_t num_write_qual = fwrite ( bufQual[i].data(), sizeof( uchar ), charsBuffered, outputFilesQual[i] );
            checkIfEqual( num_write, num_write_qual );
        }
    }
    checkIfEqual( num_write, charsBuffered );


    // closing all the output file streams
    for ( SequenceLength i = 0; i < cycleNum_; i++ )
    {
        fclose( outputFiles_[i] );
        if ( processQualities_ )
        {
            fclose( outputFilesQual[i] );
        }
    }


    lengthTexts = myTotalLen;
    lengthRead = maxSeqLen;
    kseq_destroy(seq);
    gzclose(fp);

    std::cout << "Number of sequences reading/writing: " << nSeq << "\n";
    std::cout << "Number of characters reading/writing: " << lengthTexts << "\n";

    return true;
}

/* Orders sequences by reverse lexicographic order, writes permutation of strings in file and reorders cycfiles */
bool TransposeFasta::computeRLO(const string &input, const string &output, const string &RLOsupport) {

    vector<group*> all_groups_even;
    vector<group*> all_groups_odd;

    //keeps the association between [position_of_sequence_in_file] = new_assigned_position
    keepOrderRLO = new unsigned int[nTotalSeq];
    char* cycFileContent = new char[nTotalSeq];
    //keeps the association between [position i in file] = sequence assigned at position i
    unsigned int* permutation = new unsigned int[nTotalSeq];
    bool firstTime = true;

    //TODO first iteration can be improved
    group *firstGroup = new group;
    firstGroup->base=0; firstGroup->bound=nTotalSeq-1;
    if( (cycleNum_-1) %2 == 0 )
        all_groups_even.push_back(firstGroup);
    else
        all_groups_odd.push_back(firstGroup);

    for (long long i = cycleNum_-1; i >= 0; --i) {

        vector<group*> *next_groups;
        vector<group*> *current_groups;
        if(i%2 == 0){
            next_groups = &all_groups_odd;
            current_groups = &all_groups_even;
        }else{
            next_groups = &all_groups_even;
            current_groups = &all_groups_odd;
        }

        Filename fn(output, i, "");
        outputFiles_[i] = fopen(fn, "r");
        fread(cycFileContent, sizeof(uchar), nTotalSeq, outputFiles_[i]);
        char ch;
        unsigned int position = 0;

        if(firstTime){
            //If first time read all characters and add to a unique group
            for(unsigned int h = 0; h < nTotalSeq; h++){
                charElem newChar; newChar.ch = cycFileContent[h]; newChar.oldPosition = h;
                firstGroup->newCharacters.push_back(newChar);
            }
        }else{
            //Otherwise for all groups insert its characters from cycfile
            for( group *g : *current_groups){
                //Inside group interval
                for(unsigned int b = g->base; b<= g->bound; b++) {
                    charElem newChar;
                    newChar.ch = cycFileContent[permutation[b]];
                    newChar.oldPosition = permutation[b];
                    g->newCharacters.push_back(newChar);
                }
            }
        }
        //All characters have been read and put in the right group
        //now sorting inside groups
        Filename tmp(output, "_tmp");
        FILE* fp = fopen( tmp, "w" ); //TODO check all function result
        if(fp == NULL){
            printf("Error: can't open file %s%s\n", output.c_str(), "_tmp");
            exit(EXIT_FAILURE);
        }
        unsigned int filePosition = 0;
        unsigned int realPosition = 0;

        for(unsigned int j=0; j<current_groups->size(); j++){

            group *g = (*current_groups)[j];
            //If there is a gap between two groups then print already ordered symbols
            while(filePosition != g->base){
                char chToWrite = cycFileContent[permutation[filePosition]];
                fprintf(fp, "%c", chToWrite);
                filePosition++;
            }

            sort(g->newCharacters.begin(), g->newCharacters.end(), compareCharElem);

            char prev = DUMMY_CHAR;
            char next; unsigned int howMany=0; unsigned int base = g->base;

            for(unsigned int h=0; h<g->newCharacters.size(); h++){

                //update position vector
                keepOrderRLO[g->newCharacters[h].oldPosition] = g->base+h;
                //update permutation
                permutation[g->base+h] = g->newCharacters[h].oldPosition;

                next = g->newCharacters[h].ch;
                fprintf(fp, "%c", next); //char is ordered, print it
                filePosition++;
                if((prev != next)){ //character has changed
                    //create group only if there is more than 1 element
                    if( (howMany > 1) && (prev != DUMMY_CHAR)) {
                        group *newGroup = new group;
                        newGroup->base = base;
                        newGroup->bound = base + howMany - 1;
                        next_groups->push_back(newGroup);
                    }
                    base = g->base+h;
                    howMany = 0;
                }
                howMany++;
                prev = next;
            }
            //create remaining
            if( (howMany > 1) && (prev != DUMMY_CHAR)){ //create group only if there is more than 1 element
                group *newGroup = new group;
                newGroup->base = base;
                newGroup->bound = g->bound;
                next_groups->push_back(newGroup);
            }
        }

        //write last char
        //print last char that are outside of groups
        while(filePosition < nTotalSeq){
            char chToWrite = cycFileContent[permutation[filePosition]];
            fprintf(fp, "%c", chToWrite );
            filePosition++;
        }

        //remove old cyc file
        remove(fn);
        //rename temporary file to new cyc file
        rename(tmp, fn);
        fclose(fp);

        //clean useless groups
        for(unsigned int c=0; c<current_groups->size(); c++){
            group *p = (*current_groups)[c];
            p->newCharacters.clear();
            free(p);
        } current_groups->clear();
        //no more first time
        firstTime = false;
    }

    //clean memory
    for(unsigned int c=0; c<all_groups_odd.size(); c++){
        group *p = all_groups_odd[c];
        p->newCharacters.clear();
        free(p);
    }
    for(unsigned int c=0; c<all_groups_even.size(); c++){
        group *p = all_groups_even[c];
        p->newCharacters.clear();
        free(p);
    }

    //Creating info file
    std::stringstream infoFileName;
    std::string cutNameInput = (input.substr(input.find_last_of("/\\") + 1));
    std::string::size_type const p(cutNameInput.find_last_of('.'));
    std::string file_without_extension = cutNameInput.substr(0, p);

    infoFileName << TEMP_DIR << "/" << file_without_extension.c_str() << "_RLO_permutation.txt";
    FILE* outputInfo = fopen( infoFileName.str().c_str(),"w" );
    if(!outputInfo) {
        printf("Error while opening file\n");
        exit(EXIT_FAILURE);
    }

    for(unsigned int i=0; i<nTotalSeq;i++) {
        fprintf(outputInfo, "%d\n", permutation[i]);
        fflush(outputInfo);
    }
    delete keepOrderRLO;
    delete permutation;
    fclose(outputInfo);
    return true;
}

/* ================== Not modified =======================*/
bool TransposeFasta::inputCycFile( const string &cycPrefix )
{
    //TO DO
    //The distribution of characters is useful
    //for alpha[256] -->Corresponding between the alphabet, the piles and tableOcc
    //and to know sizeAlpha

    //1) Alphabet
    //We supposed that the symbols in the input file are the following
    freq[int( terminatorChar )] = 1;
    freq[int( 'A' )] = 1;
    freq[int( 'C' )] = 1;
    freq[int( 'G' )] = 1;
    freq[int( 'N' )] = 1;
    freq[int( 'T' )] = 1;
    //GIOVANNA: ADDED THE SYMBOL Z IN THE ALPHABET, SO sizeAlpha = alphabetSize
#ifdef USE_EXTRA_CHARACTER_Z
    freq[int( 'Z' )] = 1;
#endif

    //2) Number of sequences
    string cyc1Filename = cycPrefix + "1";
    FILE *f = fopen( cyc1Filename.c_str(), "rb" );
    if ( !f )
    {
        cerr << "ERROR: Cycle file " << cyc1Filename << " not found!" << endl;
        exit( -1 );
    }
    fseek( f, 0, SEEK_END );
    nSeq = ftell( f );
    if (nSeq != ftell( f))
    {
        Logger::error() << "Error: Too many sequences. This version of BEETL was compiled for a maximum of " << maxSequenceNumber << " sequences, but this input has " << ftell(f) << " sequences. You can increase this limit by changing the type definition of 'SequenceNumber' in Types.hh and recompiling BEETL." << endl;
        exit( -1 );
    }
    fclose( f );

    //3) Length of the longest sequence
    for ( lengthRead = 1; ; ++lengthRead )
    {
        Filename cycFilename( cycPrefix, lengthRead, "" );
        FILE *f = fopen( cycFilename, "rb" );
        if ( f )
            fclose( f );
        else
            break;
    }

    //4) qualities detection
    string qual1Filename = cycPrefix + "qual.1";
    f = fopen( qual1Filename.c_str(), "rb" );
    if ( f )
    {
        processQualities_ = true;
        fclose( f );
    }
    else
        processQualities_ = false;

    //5) Total Length
    lengthTexts = lengthRead * nSeq;

    // Report
    Logger_if( LOG_SHOW_IF_VERBOSE )
    {
        Logger::out() << "****processing qualities: " << processQualities_ << "\n";
        Logger::out() << "****number of sequences: " << nSeq << "\n";
        Logger::out() << "****max length of each sequence: " << lengthRead << "\n";
        Logger::out() << "****lengthTot: " << lengthTexts << "\n";
    }

    return 1;
}
bool TransposeFasta::convertFromCycFileToFastaOrFastq( const string &fileInputPrefix, const string &fileOutput, bool generatedFilesAreTemporary, SequenceExtractor *sequenceExtractor )
{
    bool outputIsFastq = hasSuffix( fileOutput, ".fastq" );
    vector <FILE *> inFilesCyc;
    vector <FILE *> inFilesCycQual;
    //Open all cyc files
    for ( int i = 0; ; ++i )
    {
        Filename fn( fileInputPrefix, i, "" );
        FILE *f = fopen( fn, "rb" );
        if ( !f ) break;
        inFilesCyc.push_back( f );

        if ( outputIsFastq )
        {
            Filename fnQual( fileInputPrefix, i, ".qual" );
            inFilesCycQual.push_back( fopen( fnQual, "rb" ) );
            if ( inFilesCycQual[i] == NULL )
            {
                std::cerr << "TransposeFasta: could not open file " << fnQual << std::endl;
                exit ( EXIT_FAILURE );
            }
        }
    }
    if ( inFilesCyc.empty() )
    {
        std::cerr << "TransposeFasta: could not open file " << fileInputPrefix << "0" << std::endl;
        exit ( EXIT_FAILURE );
    }
    SequenceLength lengthRead = inFilesCyc.size();
    fseek( inFilesCyc[0], 0, SEEK_END );
    SequenceNumber nSeq = ftell( inFilesCyc[0] );
    fseek( inFilesCyc[0], 0, SEEK_SET );

    ofstream outFile ( fileOutput.c_str() );
    if ( outFile.is_open() == false )
    {
        std::cerr << "Error opening \"" << fileOutput << "\" file" << std::endl;
        exit ( 1 );
    }

    //I must read a char for each sequence. The chars at the position i corresponds to the chars of the sequence i.
    char symbol;
    string sequence = "";
    // buf to accelerate SequenceExtractor usage
    const int SEQ_EXTRACTION_BUF_SIZE = 1024;
    char seqExtractionBuf[SEQ_EXTRACTION_BUF_SIZE];
    int seqCountToSkip = 0;
    for ( SequenceNumber j = 0; j < nSeq; j++ )
    {
        bool extractThisSeq = !sequenceExtractor || sequenceExtractor->doWeExtractNextSequence();

        if ( !extractThisSeq )
        {
            ++seqCountToSkip;
            continue;
        }
        else
        {
            while ( seqCountToSkip > 0 )
            {
                size_t skip = min( seqCountToSkip, SEQ_EXTRACTION_BUF_SIZE );
                for ( SequenceLength i = 0; i < lengthRead; i++ )
                {
                    assert( fread ( seqExtractionBuf, sizeof( char ), skip, inFilesCyc[i] ) == skip );
                    if ( outputIsFastq && inFilesCycQual.size() >= lengthRead )
                        assert( fread ( seqExtractionBuf, sizeof( char ), skip, inFilesCycQual[i] ) == skip );
                }
                seqCountToSkip -= skip;
            }
        }

        if ( outputIsFastq )
            outFile << "@Read"  << j << std::endl;
        else
            outFile << "> Read "  << j << std::endl;
        for ( SequenceLength i = 0; i < lengthRead; i++ )
        {
            assert( fread ( &symbol, sizeof( char ), 1, inFilesCyc[i] ) == 1 );
            sequence.append ( 1, symbol );
        }
        outFile << sequence << std::endl;
        Logger_if( LOG_FOR_DEBUGGING ) Logger::out() << sequence << std::endl;
        sequence.clear();

        if ( outputIsFastq )
        {
            outFile << "+" << std::endl;
            if ( outputIsFastq && inFilesCycQual.size() >= lengthRead )
            {
                for ( SequenceLength i = 0; i < lengthRead; i++ )
                {
                    assert( fread ( &symbol, sizeof( char ), 1, inFilesCycQual[i] ) == 1 );
                    sequence.append ( 1, symbol );
                }
                outFile << sequence << std::endl;
                sequence.clear();
            }
            else
                outFile << "<qualities not available>" << std::endl;
        }
    }


    outFile.close();


    //Close all cyc files
    for ( SequenceLength i = 0; i < lengthRead; i++ )
    {
        fclose( inFilesCyc[i] );
    }

    return 1;
}


/* ================== Dismissed methods =================*/
//Sort by len
/*void TransposeFasta::sortBylen( const string &input, const string &output){ //TODO add qualities

    unsigned int contaDollari = 0;
    unsigned int contaCancelletti = 0;
    unsigned int* daSeqAPosizione = new unsigned int [nTotalSeq];
    vector <unsigned int> nextIterazione;
    vector <unsigned int> thisIterazione;
    char* symbToPrint = new char[nTotalSeq];
    char ch; //char to read
    unsigned int i = 0;
    unsigned int contaLettere = 0;
    unsigned int posizione = 0;
    //First iteration

    //Open cyc file 0
    Filename fn(output, 0, "");
    outputFiles_[i] = fopen(fn, "r");
    if(outputFiles_[i] == NULL){
        cerr<<"Error: can't open file cyc.0"<<endl;
        exit(EXIT_FAILURE);
    }

    //Open tmp cyc file
    Filename tmp(output, "_tmp");
    FILE* fp = fopen( tmp, "w" );
    if(fp == NULL){
        cerr<<"Error: can't open file"<<output<<"_tmp";
        exit(EXIT_FAILURE);
    }

    //I read one char at time
    unsigned int sequenceIndex = 0;
    while ((ch = fgetc(outputFiles_[i])) != EOF){

        if(ch == DUMMY_CHAR)
            contaCancelletti++;

        if(ch == TERMINATE_CHAR) {
            contaDollari++;
            nextIterazione.push_back(sequenceIndex);
        }

        //I establish ordering
        if( (ch != DUMMY_CHAR) && (ch != TERMINATE_CHAR) ) {
            //it's a letter A,C,G,N,T I have to store its pos
            thisIterazione.push_back(sequenceIndex);
            //scrivo char nel file
            fprintf(fp, "%c", ch);
        }

        sequenceIndex++;
    }

    //write all $
    for(unsigned int dollar = 0; dollar < contaDollari; ++dollar){
        unsigned int printed = fprintf(fp, "$");
    }

    //write all #
    for(unsigned int canc = 0; canc < contaCancelletti; ++canc){
        unsigned int printed = fprintf(fp, "#");
    }

    remove(fn);
    //rename temporary file to new cyc file
    rename(tmp, fn);
    fclose(fp);

    //merge thisIterazione + nextIterazione, clean nextIterazione
    thisIterazione.insert(thisIterazione.end(), nextIterazione.begin(), nextIterazione.end()),
    nextIterazione.clear();


    //From second to cycleNum-1 iteration

    for (i = 1; i < cycleNum_; ++i) {
        contaLettere = 0;
        contaDollari = 0;
        contaCancelletti = 0;

        //Open cyc file
        Filename fn(output, i, "");
        outputFiles_[i] = fopen(fn, "r");
        if(outputFiles_[i] == NULL){
            cerr<<"Error: can't open file"<<output<<endl;
            exit(EXIT_FAILURE);
        }

        //Open tmp cyc file
        Filename tmp(output, "_tmp");
        FILE* fp = fopen( tmp, "w" );
        if(fp == NULL){
            cerr<<"Error: can't open file"<<output<<"_tmp";
            exit(EXIT_FAILURE);
        }

        //I read one char at time
        unsigned int sequenceIndex = 0;
        while ((ch = fgetc(outputFiles_[i])) != EOF){

            if(ch == DUMMY_CHAR)
                contaCancelletti++;

            if(ch == TERMINATE_CHAR) {
                contaDollari++;
                nextIterazione.push_back(sequenceIndex);
            }

            //I establish ordering
            if( (ch != DUMMY_CHAR) && (ch != TERMINATE_CHAR) ) {
                symbToPrint[ thisIterazione[contaLettere] ] = ch;
                printf("Added %c\n",symbToPrint[thisIterazione[contaLettere]]);
                contaLettere++;
            }
            sequenceIndex++;
        }

        //print all symbToPrint
        size_t num_write_len = fwrite(symbToPrint, sizeof(char), thisIterazione.size(), fp);

        //write all $
        for(unsigned int dollar = 0; dollar < contaDollari; ++dollar){
            unsigned int printed = fprintf(fp, "$");
        }

        //write all #
        for(unsigned int canc = 0; canc < contaCancelletti; ++canc){
            unsigned int printed = fprintf(fp, "#");
        }

        remove(fn);
        //rename temporary file to new cyc file
        rename(tmp, fn);
        fclose(fp);

        //merge thisIterazione + nextIterazione, clean nextIterazione
        thisIterazione.insert(thisIterazione.end(), nextIterazione.begin(), nextIterazione.end()),
        nextIterazione.clear();

    }

    delete[] symbToPrint;

}
*/

//Sort by len 2
/* Creates cyc files with strings ordered by len and filling difference in size with DUMMY_CHAR, prints permutation of strings in file
bool TransposeFasta::convertByLen(const string &input, const string &output, bool generatedFilesAreTemporary ){

    if(differentLenDetected_ == false) {
        Logger::out() << "Asked to order sequences by length but all sequences are the same size. Sorting by length canceled."<<endl;
        return convert(input, output, generatedFilesAreTemporary);
    }

    vector<vector<uchar> > bufQual;
    vector<FILE *> outputFilesQual;

    if ( processQualities_ )
    {
        bufQual.resize( cycleNum_, vector<uchar>( BUFFERSIZE ) );
        outputFilesQual.resize( cycleNum_ );
    }

    //TO DO
    lengthRead = cycleNum_;
    //The distribution of characters is useful
    //for alpha[256] -->Corresponding between the alphabet, the piles and tableOcc
    //and to know sizeAlpha
    //We supposed that the symbols in the input file are the following
    freq[int( terminatorChar )] = 1;
    freq[int( 'A' )] = 1;
    freq[int( 'C' )] = 1;
    freq[int( 'G' )] = 1;
    freq[int( 'N' )] = 1;
    freq[int( 'T' )] = 1;
    //GIOVANNA: ADDED THE SYMBOL Z IN THE ALPHABET, SO sizeAlpha = alphabetSize
#ifdef USE_EXTRA_CHARACTER_Z
    freq[int( 'Z' )] = 1;
#endif

    // create output files
    for ( SequenceLength i = 0; i < cycleNum_; i++ )
    {
        Filename fn( output, i, "" );
        outputFiles_[i] = fopen( fn, "w" );
        if ( outputFiles_[i] == NULL )
        {
            cerr << "Error: couldn't open output file " << fn << endl;
            if ( i > 0 )
            {
                cerr << "  You may have reached the maximum number of opened files (see `ulimit -n`) or the maximum number of files allowed in one directory, as we create one file per cycle (and a second one if qualities are present)" << endl;
                exit ( -1 );
            }
        }
        if ( generatedFilesAreTemporary )
            TemporaryFilesManager::get().addFilename( fn );
        if ( processQualities_ )
        {
            Filename fnQual( output + "qual.", i, "" );
            outputFilesQual[i] = fopen( fnQual, "w" );
            if ( outputFilesQual[i] == NULL )
            {
                cerr << "Error: couldn't open output file " << fnQual << endl;
                if ( i > 0 )
                {
                    cerr << "  You may have reached the maximum number of opened files (see `ulimit -n`) or the maximum number of files allowed in one directory, as we create one file per cycle (and a second one if qualities are present)" << endl;
                    exit ( -1 );
                }
            }
            if ( generatedFilesAreTemporary )
                TemporaryFilesManager::get().addFilename( fnQual );
        }
    }

    //calculating how many times must read input file
    double timesReadInput = (double)nTotalSeq/(double)BUFFERSIZE;
    if(timesReadInput - (unsigned)timesReadInput > 0) timesReadInput = (unsigned)timesReadInput + 1;
    //printf("\n Sequences are %ld and Times to read found is %f\n",nTotalSeq, (double)timesReadInput);
    int counterTimes = 0;
    nSeq = 0;
    int base = 0;
    int bound = BUFFERSIZE - 1;
    lengthTexts = 0;
    unsigned int num_write=0;
    unsigned int charsBuffered = 0;
    int maxPosition = 0;
    int insertedSequences = 0;

    // looping through the input file, add the characters to the buffer, print buffer when it's full
    while( (counterTimes < timesReadInput) && (insertedSequences < nTotalSeq) ) {

        num_write = 0;
        nSeq=0;
        gzFile fp;
        kseq_t *seq;
        int len = 0;
        fp = gzopen(input.c_str(), "r");
        seq = kseq_init(fp);
        maxPosition = 0;

        //TODO uso buffer o no? Se ne scrivo una alla volta?
        while ((len = kseq_read(seq)) >= 0 && (insertedSequences < nTotalSeq)) {

            int position = keepOrder[nSeq];
            //printf("\n File to open is %s Position found is %d\n", input.c_str(), position);
            if (charsBuffered == BUFFERSIZE) {
                // write buffers to the files, clear buffers
#pragma omp parallel for num_threads(4)
                for (int i = 0; i < cycleNum_; i++) {
                    size_t num_write_bases = fwrite(buf_[i].data(), sizeof(char), charsBuffered, outputFiles_[i]);
                    checkIfEqual(num_write_bases,
                                 charsBuffered); // we should always read/write the same number of characters
                    if (processQualities_) {
                        size_t num_write_qual = fwrite(bufQual[i].data(), sizeof(char), charsBuffered,
                                                       outputFilesQual[i]);
                        checkIfEqual(num_write_bases, num_write_qual);
                    }
                }
                //lengthTexts += (num_write * cycleNum_);
                charsBuffered = 0;
                maxPosition = 0;
                base += BUFFERSIZE;
                bound += BUFFERSIZE;

                /* debug buffer print
                printf("\n--------------------------\n");
                for(int i = 0; i < cycleNum_; i++) {
                  for (int j = 0; j < BUFFERSIZE; j++)
                      printf("%c ", buf_[i][j]);
                  printf("\n");
                }
                printf("\n--------------------------\n");


                //reset buffer
                for(int i = 0; i < cycleNum_; i++)
                    for(int j = 0; j < BUFFERSIZE; j++)
                        buf_[i][j] = DUMMY_CHAR;
            }

            //check if sequence has to be ignored for now
            if ((position >= base) && (position <= bound)) {
                int bufferPosition = position % BUFFERSIZE;
                if(bufferPosition > maxPosition)
                    maxPosition = bufferPosition;

#if (ALIGN == 0)
                //Align right side to compute preprocessing RLO
                int index = cycleNum_-1;
                for ( int i = len-1; i >= 0; --i ){
                    buf_[index][bufferPosition] = seq->seq.s[i];

                    if (processQualities_ ){
                        bufQual[index][bufferPosition] = seq->qual.s[i];
                    }
                    index--;
                }
                if(index >= 0)
                    buf_[index][bufferPosition] = TERMINATE_CHAR;
#else
                //else align left side
                SequenceLength i;
                for ( i = 0; i < len; ++i ) {
                    buf_[i][bufferPosition] = seq->seq.s[i];

                    if (processQualities_) {
                        bufQual[i][bufferPosition] = seq->qual.s[i];
                    }
                }
                 if(i <= cycleNum_-1)
                    buf_[i][bufferPosition] = TERMINATE_CHAR;
#endif

                // increase the counter of chars buffered only if the sequence was inserted
                charsBuffered++;
                // increase the number of inserted sequences
                insertedSequences++;
                lengthTexts+=len;
            }
            nSeq++;
        }
        counterTimes++;
    }

    // write the rest
    for (SequenceLength i = 0; i < cycleNum_; i++) {
        num_write = fwrite(buf_[i].data(), sizeof(uchar), maxPosition+1, outputFiles_[i]);
        //lengthTexts += num_write;
        if (processQualities_) {
            size_t num_write_qual = fwrite(bufQual[i].data(), sizeof(uchar), maxPosition+1, outputFilesQual[i]);
            checkIfEqual(num_write, num_write_qual);
        }
    }
    checkIfEqual(num_write, maxPosition+1);

    // closing all the output file streams
    for ( SequenceLength i = 0; i < cycleNum_; i++ )
    {
        fclose( outputFiles_[i] );
        if ( processQualities_ )
        {
            fclose( outputFilesQual[i] );
        }
    }

    std::cout << "Number of sequences reading/writing: " << nSeq << "\n";
    std::cout << "Number of characters reading/writing: " << lengthTexts << "\n";

    return true;
}
*/

//Find RLo group
/*int findGroup(int realPosition, vector<group*> &RLOgroups){
//TODO use iterator and improve, if position > bound -> then no group
    int groupID=-1;
    int iter=0;

    for(group *g : RLOgroups){

        if (iter > g->bound) //se ho superato il bound del gruppo posso uscire, i gruppi sono ordinati per costruzione
            break;

        if ( (realPosition >= g->base) && (realPosition <= g->bound) ) {
            groupID = iter;
            break;
        }
        iter++;
    }
    return groupID;
}*/

