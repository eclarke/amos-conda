#include <iostream>

#include "Align_Utils.hh"

namespace Output 
{

static const int32_t ERROR_FACTOR = 10;
static const uint32_t MAX_SIZE = 1 << 31;

bool hasValidOverlap(
        const std::string sequenceA,
        const std::string sequenceB,
        uint32_t offset,
        uint32_t startA,
        uint32_t endA,
        uint32_t startB,
        uint32_t endB,
        bool     &contained,
        Alignment_t &ali) {
   const std::string* first = NULL;
   const std::string* second = NULL;
   contained = false;
   std::string overlap;

   int startFirst, endFirst, startSecond, endSecond;
   int maxPos = 0;
   if (startA < startB) {
      first = &sequenceA;
      startFirst = startA;
      endFirst = endA;
      second = &sequenceB;
      startSecond = startB;
      endSecond = endB;
      maxPos = endA-startA;
   } else {
      first = &sequenceB;                
      startFirst = startB;
      endFirst = endB;
      second = &sequenceA;                
      startSecond = startA;
      endSecond = endA;
      maxPos = endB - startB;
   }
   uint32_t ovlStart = (startA > startB ? startA : startB);
   uint32_t ovlEnd = (endA < endB ? endA : endB);
   if ((int)ovlEnd - (int)ovlStart > 0) {
      if (first->length() < second->length()) {
         if ((int)startFirst - (int)startSecond + OVL_WIGGLE >= 0 && (int)endSecond - (int)endFirst + OVL_WIGGLE >= 0) {
            contained = true;
            overlap = second->substr(ovlStart - Min(startA, startB), ovlEnd - ovlStart);
         }
      } else {
         if ((int)startSecond - (int)startFirst + OVL_WIGGLE >= 0 && (int)endFirst - (int)endSecond + OVL_WIGGLE >= 0) {
            contained = true;
            overlap = first->substr(ovlStart - Min(startA, startB), ovlEnd - ovlStart);
         }
      }
   }

   bool aligns = false;
   int32_t numErrors = (int32_t) (roundf((Output::ERROR_RATE)*((float)ovlEnd-ovlStart)));
   int32_t minOvl = (int)roundf(((float)ovlEnd-ovlStart)/2);
   if (minOvl < 0 && ovlStart >= ovlEnd) { minOvl = Output::MIN_OVERLAP; }
   if (numErrors < 0 && ovlStart >= ovlEnd) { numErrors = (int32_t)(roundf((Output::ERROR_RATE)*2*Output::OVL_WIGGLE)); }
   int32_t gapWiggle = (ovlStart >= ovlEnd ? Output::OVL_WIGGLE*2 : Output::OVL_WIGGLE);

   if (contained) {
      int count = 0;
      for (int i = 0; i < overlap.length(); i++) {
         if (overlap[i] == SCAFFOLD_GAP_CHAR) { count++; }
      }
      if ((double)count / overlap.length() > 0.90) {
         // no need to do alignment just replace the region, return alignment info
         ali.a_lo = ovlStart - Min(startA, startB); 
         ali.a_hi = ovlEnd - ovlStart + ali.a_lo;
         ali.b_lo = ali.a_lo; 
         ali.b_hi = ali.a_hi; 
         ali.errors = 0;
         return true;
      }
   }

   while (aligns == 0 && minOvl >= Output::MIN_OVERLAP) {
       int32_t loPos = Max(0, (int32_t)Min(ovlStart, ovlEnd) - (int32_t)offset - gapWiggle);
       int32_t hiPos = Min(maxPos, loPos + ERROR_FACTOR*Output::OVL_WIGGLE);
       int32_t wiggleArea = hiPos - loPos;
       // the overlap match code cannot handle more than a 2^31 total array, lower error rate to fit
       double d = (double)(1 + numErrors) * (double)(numErrors + 1 + wiggleArea);
       if (d > MAX_SIZE) {
           numErrors = (int)floor(sqrt(MAX_SIZE) - ((double)wiggleArea/2));
       }
      aligns = Overlap_Match_VS(second->c_str(), second->length(), first->c_str(), first->length(),
                           loPos,
                           hiPos,//maxPos, //Min(maxPos, (int)Min(ovlStart, ovlEnd) - (int32_t)offset + gapWiggle),
                           minOvl/*Output::MIN_OVERLAP*/, numErrors, ali);

       // check if this overlap is acceptable
       int32_t ovlSize = Min(ali.a_hi - ali.a_lo, ali.b_hi - ali.b_lo);
       if (ovlSize < Output::MIN_OVERLAP || ali.errors > Output::ERROR_RATE*(float)(ovlSize)) {
          aligns = false;
       }
       if (contained && (ali.a_hi-ali.a_lo != second->length() && ali.b_hi-ali.b_lo != first->length())) {
          aligns = false; 
       }
       if (minOvl > Output::MIN_OVERLAP) {
          minOvl = minOvl / 2;
          if (minOvl < Output::MIN_OVERLAP) { minOvl = Output::MIN_OVERLAP; }
       } else {
          break;
       }
    }

    return aligns;
}
}
