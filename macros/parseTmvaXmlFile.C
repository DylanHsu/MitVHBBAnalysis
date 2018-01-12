#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

string getXmlProperty( string line, string property ) {
  string token=property+"=\"";
  string value="";
  size_t pos1, pos2;
  pos1 = line.find(token); if(pos1==string::npos) return value;
  pos2 = line.find("\"", pos1+token.length()); if(pos2==string::npos) return value;
  value=line.substr(pos1+token.length(), pos2-pos1-token.length());
  return value;
}

using namespace std;
// Returns vector of (Label, Expression) pairs
vector<pair<string, string> > parseTmvaXmlFile(
  string xmlFile
) {
  ifstream ifs; ifs.open(xmlFile.c_str());
  assert(ifs.is_open());
  string line;
  bool startOfList=false;
  unsigned nVar=0;
  vector<pair<string, string> > varLabelExprs;
  while(getline(ifs,line)) {
    size_t pos1,pos2; string token;
    
    if(!startOfList) {
      // <Variables NVar="34">
      string nVarStr = getXmlProperty(line, "NVar");
      if(nVarStr=="") continue;
      nVar = atoi( nVarStr.c_str());
      assert(nVar!=0); // "" => 0
      startOfList=true;
    }
    if(line.find("</Variables>") != string::npos) break;
    // <Variable VarIndex="8" Expression="mT" Label="mT" Title="W boson m_{T}" Unit="GeV" Internal="mT" Type="F" Min="2.74883292e-04" Max="4.35341675e+02"/>
    string varIndexStr = getXmlProperty(line, "VarIndex");
    string varExpr = getXmlProperty(line, "Expression");
    string varLabel = getXmlProperty(line, "Label");
    if(varIndexStr=="" || varIndexStr=="" || varLabel=="") continue;
    varLabelExprs.emplace_back(varLabel, varExpr);
    
    printf("Parsed variable #%d: label \"%s\" => expression \"%s\"\n", atoi(varIndexStr.c_str()), varLabel.c_str(), varExpr.c_str());
    if(varLabelExprs.size()>=nVar || line.find("</Variables>") != string::npos)  break;
  }
  ifs.close();
  return varLabelExprs;
}
