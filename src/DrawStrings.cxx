#include "DrawStrings.h"


/** 
 * Constructs a new cut, where the result is multiplied by "weight"
 * 
 * @param cut is the cut whose result should be multiplied by weight (default is empty cut)
 * 
 * @return the new weighted cut
 */
const TCut Acclaim::ThermalTree::weight(const TCut cut){

  if(strlen(cut.GetTitle()) > 0){
    return TCut(TString::Format("weight*(%s)", cut.GetTitle()).Data());
  }
  else{
    return TCut("weight");
  }
}
