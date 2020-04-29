/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_compat_2_h
#define _aspect_compat_2_h

#include <aspect/global.h>
#include <aspect/utilities.h>

/**
 * Add a function to Utilities that computes the distance
 * of a point to a line segment.
 */
namespace aspect
{
  namespace Utilities
  { 
    /**  
     * Given a 2d point and a list of two points that define a line, compute the smallest
     * distance of the point to the line segment. Note that when the point's perpendicular
     * base does not lie on the line segment, the smallest distance to the segment's end
     * points is calculated.
     */
    template <int dim> 
    double
    distance_to_line(const std::vector<dealii::Point<2> > &point_list,
                     const dealii::Point<2> &point);
  }

}

namespace aspect
{
    template <int dim>
    inline 
    double
    distance_to_line(const std::vector<dealii::Point<2> > &point_list,
                     const dealii::Point<2> &point)
    {    

      /**  
       * This code is based on http://geomalgorithms.com/a02-_lines.html#Distance-to-Infinite-Line,
       * and therefore requires the following copyright notice:
       *
       * Copyright 2000 softSurfer, 2012 Dan Sunday
       * This code may be freely used and modified for any purpose
       * providing that this copyright notice is included with it.
       * SoftSurfer makes no warranty for this code, and cannot be held
       * liable for any real or imagined damage resulting from its use.
       * Users of this code must verify correctness for their application.
       *
       */

      const unsigned int n_poly_points = point_list.size();
      AssertThrow(n_poly_points == 2, ExcMessage("A list of points for a line segment should consist of 2 points."));

      // Create vector along the polygon line segment P0 to P1
      Tensor<1,2> vector_segment = point_list[1] - point_list[0];
      // Create vector from point P to the second segment point
      Tensor<1,2> vector_point_segment = point - point_list[0];

      // Compute dot products to get angles
      const double c1 = vector_point_segment * vector_segment;

      // Point P's perpendicular base line lies outside segment, before P0.
      // Return distance between points P and P0.
      if (c1 <= 0.0) 
        return (Tensor<1,2> (point_list[0] - point)).norm();

      const double c2 = vector_segment * vector_segment;

      // Point P's perpendicular base line lies outside segment, after P1.
      // Return distance between points P and P1.
      if (c2 <= c1)
        return (Tensor<1,2> (point_list[1] - point)).norm();

      // Point P's perpendicular base line lies on the line segment.
      // Return distance between point P and the base point.
      const Point<2> point_on_segment = point_list[0] + (c1/c2) * vector_segment;
      return (Tensor<1,2> (point - point_on_segment)).norm();

    }
}


#endif
