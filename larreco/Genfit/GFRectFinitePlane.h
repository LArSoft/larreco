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

/**
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */

/** @addtogroup genfit
 * @{
 */

#ifndef GFRECTFINITEPLANE_H
#define GFRECTFINITEPLANE_H

#include "larreco/Genfit/GFAbsFinitePlane.h"
#include <iostream>
#include <stdexcept> // std::logic_error

/** @brief Concrete implementation of finitie detector plane for rectangles.
 */

namespace genf {

  class GFRectFinitePlane : public GFAbsFinitePlane {
  public:
    //override inActive & Print methods
    bool inActive(const double& u, const double& v) const;
    void Print(std::ostream& out = std::cout) const;

    //! give dimensions of finite rectangle: u1,u2,v1,v2
    GFRectFinitePlane(const double&, const double&, const double&, const double&);
    GFRectFinitePlane();

    virtual ~GFRectFinitePlane();

    GFAbsFinitePlane* clone() const { return new GFRectFinitePlane(*this); }

  private:
    double fUmin, fUmax, fVmin, fVmax;

    virtual void Print(Option_t*) const
    {
      throw std::logic_error(std::string(__func__) + "::Print(Option_t*) not available");
    }
    // public:
    //ClassDef(GFRectFinitePlane,1)
  };

} // namespace genf

#endif

/** @} */
