/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
/** @addtogroup genfit
 * @{ */

#ifndef GFEXCEPTION_H
#define GFEXCEPTION_H

#include <exception>
#include <iomanip> // std::setw()
#include <ios>     // std::ios::fmtflags
#include <sstream>
#include <string>
#include <vector>

#include "RtypesCore.h"
#include "TMatrixT.h"
#include "TVector3.h"

/** @brief Exception class for error handling in GENFIT (provides storage for diagnostic information)
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * This is the class that is used for all error handling in GENFIT.
 * It is a utility class that allows to store numbers and matrices together
 * with an error string. The exception class can then be thrown when an error
 * is detected and the C++ exception handling facilities can be used to
 * catch and process the exception.
 */
class GFException : public std::exception {
private:
  static bool fQuiet;

  std::string fExcString;
  int fLine;
  std::string fFile;

  std::string fNumbersLabel;
  std::string fMatricesLabel;
  std::vector<double> fNumbers;
  std::vector<TMatrixT<Double_t>> fMatrices;

  bool fFatal;

public:
  /** @brief Initializing constructor
   *
   * @param what error message
   * @param line line at which the exception is created. Can be set through
   * __LINE__ macro
   * @param file sorcefile in which the exception is created.
   * Can be set through __FILE__ macro
   */
  GFException(std::string, int, std::string);
  virtual ~GFException() throw();

  /** @brief set fatal flag. if this is true, the fit stops for this current track repr. */
  GFException& setFatal(bool b = true)
  {
    fFatal = b;
    return *this;
  }
  /** @brief get fatal flag. */
  bool isFatal() { return fFatal; }
  /** @brief set list of numbers with description */
  GFException& setNumbers(std::string, const std::vector<double>&);
  /** @brief set list of matrices with description */
  GFException& setMatrices(std::string, const std::vector<TMatrixT<Double_t>>&);

  /** @brief print information in the exception object */
  void info();

  //! standard error message handling for exceptions. use like "std::cerr << e.what();"
  const char* what() const noexcept override;

  std::string getExcString() { return fExcString; }

  static void quiet(bool b = true) { fQuiet = b; }
};

namespace genf {
  //------------------------------------------------------------------------------
  //@{
  /// Small utility functions which print some ROOT objects into an output stream
  template <class ROOTOBJ>
  void PrintROOTobject(std::ostream&, const ROOTOBJ&);

  template <>
  void PrintROOTobject(std::ostream&, const TVector3& v);

  template <typename T>
  void PrintROOTmatrix(std::ostream& out, const TMatrixT<T>& m);
  //@}

  /// Shortcut to write one ROOT object into a string
  template <class ROOTOBJ>
  std::string ROOTobjectToString(const ROOTOBJ& obj)
  {
    std::ostringstream sstr;
    PrintROOTobject(sstr, obj);
    return sstr.str();
  }

} // namespace genf

//------------------------------------------------------------------------------
// template definitions
//
template <class ROOTOBJ>
void genf::PrintROOTobject(std::ostream&, const ROOTOBJ& obj)
{
  obj.Print();
}

template <typename T>
void genf::PrintROOTmatrix(std::ostream& out, const TMatrixT<T>& m)
{

  constexpr std::streamsize fw = 11;
  constexpr std::streamsize ifw = 4 + (fw & 1);
  const Int_t rb = m.GetRowLwb(), cb = m.GetColLwb();

  const Int_t R = m.GetNrows(), C = m.GetNcols();
  out << R << "x" << C << " matrix is as follows";

  std::streamsize swidth = out.width(4);
  std::ios::fmtflags sflags = out.flags();
  out.unsetf(std::ios_base::floatfield); // out << std::defaultfloat;

  // header: column number
  std::string index_pad((fw - ifw) / 2, ' ');
  out << "\n" << std::string(ifw, ' ') << " |";
  for (Int_t c = 0; c < C; ++c)
    out << index_pad << std::setw(ifw) << (cb + c) << index_pad << "|";

  // dashed line
  out << "\n" << std::string((C + 1) * (fw + 1), '-');

  // content, row by row
  for (Int_t r = 0; r < R; ++r) {
    // header: row number
    out << "\n" << std::setw(ifw) << (rb + r) << " |";
    for (Int_t c = 0; c < C; ++c)
      out << std::setw(fw) << m(rb + r, cb + c) << " ";
  } // for r
  out << "\n\n";

  // restore the stream features
  out.flags(sflags);
  out.width(swidth);
} // genf::PrintROOTmatrix<TMatrixT<T>>()

#endif

/** @} */
