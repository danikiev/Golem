/******************************************************************************/
/*           GOLEM - Multiphysics of faulted geothermal reservoirs            */
/*                                                                            */
/*          Copyright (C) 2017 by Antoine B. Jacquey and Mauro Cacace         */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*    This program is free software: you can redistribute it and/or modify    */
/*    it under the terms of the GNU General Public License as published by    */
/*      the Free Software Foundation, either version 3 of the License, or     */
/*                     (at your option) any later version.                    */
/*                                                                            */
/*       This program is distributed in the hope that it will be useful,      */
/*       but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       */
/*                GNU General Public License for more details.                */
/*                                                                            */
/*      You should have received a copy of the GNU General Public License     */
/*    along with this program.  If not, see <http://www.gnu.org/licenses/>    */
/******************************************************************************/
#ifndef GOLEMMAP_H
#define GOLEMMAP_H

#include "ElementUserObject.h"

class GolemMap;

template <>
InputParameters validParams<GolemMap>();

class GolemMap : public ElementUserObject
{
public:
  GolemMap(const InputParameters & parameters);
  ~GolemMap() {}

  virtual void initialSetup() override;

  virtual void initialize() {}
  virtual void execute() {}
  virtual void finalize() {}
  virtual void threadJoin(const UserObject & y) {}

protected:
  bool mapElements(const unsigned int subdomain_id);
  bool pointsAreEqual(const Point pt0, const Point pt1, const Real tol);
  void writeFile();

  FileName _file_name;

  std::vector<unsigned int> _frac_elem_id;
  std::vector<unsigned int> _matrix_elem_id;
  const Real _tolerance;
  std::map<unsigned int, std::vector<unsigned int>> _mapped_elem_id;
};

#endif // GOLEMMAPMATRIXTOFRAC_H
