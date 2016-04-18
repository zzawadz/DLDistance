#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
using namespace Rcpp;

std::vector<unsigned int> traverse_h(const std::vector<std::vector<unsigned int> >& H, bool preferSubs)
{
  unsigned int insCntr  = 0;
  unsigned int subsCntr = 0;
  unsigned int delCntr  = 0;
  unsigned int swapCntr = 0;

  size_t i = H.size() - 1;
  size_t j = H.at(0).size() -1 ;

  while(H[i][j] > 0)
  {
    unsigned int del  = H[i-1][j];
    unsigned int ins  = H[i][j-1];
    unsigned int subs = H[i - 1][j - 1];

    if(i > 1 && j > 1)
    {
      unsigned int swap = H[i - 2][j - 2];

      unsigned int hij = H[i][j];


      if(hij - swap == 1 &&
         hij == del &&
         hij == ins &&
         hij == subs)
      {
        swapCntr++;
        i-=2;
        j-=2;
        continue;
      }
    }


    unsigned int val = std::min(std::min(subs,ins),del);

    if(preferSubs)
    {
      if(val == subs)
      {
        if(H[i - 1][j - 1] != H[i][j])
        {
          subsCntr++;
        }
        i--;
        j--;

      } else
        if(val == ins)
        {
          j--;
          insCntr++;

        } else if(val == del)
        {
          i--;
          delCntr++;
        }
    }
    else{
      if(val == ins)
      {
        j--;
        insCntr++;

      } else if(val == del)
      {
        i--;
        delCntr++;
      } else if(val == subs)
      {

        if(H[i - 1][j - 1] != H[i][j])
        {
          subsCntr++;
        }
        i--;
        j--;

      }
    }

  }

  std::vector<unsigned int> result(4);
  result[0] = insCntr;
  result[1] = delCntr;
  result[2] = subsCntr;
  result[3] = swapCntr;


  return result;
}



//' @export
// [[Rcpp::export]]
std::vector<unsigned int> dl_distance(std::string a, std::string b) {


  if(a == b) return std::vector<unsigned int>(4);

  unsigned int inf = a.size() + b.size();

  std::vector<std::vector<unsigned int> > H(a.size() + 2, std::vector<unsigned int>(b.size()+2));

  H[0][0] = inf;

  for(unsigned int i = 0; i<=a.length(); i++) {
    H[i+1][1] = i;
    H[i+1][0] = inf;
  }
  for(unsigned int j = 0; j<=b.length(); j++) {
    H[1][j+1] = j;
    H[0][j+1] = inf;
  }

  std::vector<unsigned int> DA(256);



  for(unsigned int i = 1; i <= a.length(); i++)
  {
    unsigned int DB = 0;
    for(unsigned int j = 1; j <= b.length(); j++) {
      unsigned int i1 = DA[b.at(j-1)];
      unsigned int j1 = DB;
      unsigned int d = ((a.at(i-1)  ==b.at(j-1) ) ? 0:1);
      if(d==0) DB = j;

      unsigned int subs = H[i][j]+d;
      unsigned int ins  = H[i+1][j] + 1;
      unsigned int del  = H[i][j+1] + 1;
      unsigned int swap = H[i1][j1] + (i-i1-1) + 1 + (j-j1-1);

      unsigned int val =
        std::min(std::min(std::min(subs,
                                   ins),
                                   del),
                                   swap);



      H[i+1][j+1] = val;
    }

    DA[a.at(i-1)] = i;

  }


  std::vector<unsigned int> res1 = traverse_h(H, false); // prefer deletions
  std::vector<unsigned int> res2 = traverse_h(H, true); // prefer Substitutions

  unsigned int r1 = std::accumulate(res1.begin(), res1.end(), 0);
  unsigned int r2 = std::accumulate(res2.begin(), res2.end(), 0);

  unsigned int dist = H.back().back();

  if(dist == r1) return res1;
  if(dist == r2) return res2;

  Rf_error("Cannot match D-L distance with performed operations");

  return std::vector<unsigned int>(4);
}



