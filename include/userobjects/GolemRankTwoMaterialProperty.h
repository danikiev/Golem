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
#ifndef GOLEMRANKTWOMATERIALPROPERTY_H
#define GOLEMRANKTWOMATERIALPROPERTY_H

#include "ElementIntegralUserObject.h"

class GolemRankTwoMaterialProperty;

template <>
InputParameters validParams<GolemRankTwoMaterialProperty>();

class GolemRankTwoMaterialProperty : public ElementIntegralUserObject
{
public:
  GolemRankTwoMaterialProperty(const InputParameters & parameters);
  virtual ~GolemRankTwoMaterialProperty() {}

  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual void threadJoin(const UserObject & y);

  virtual Real computeIntegral();
  virtual Real computeQpIntegral();

  Real getElementalValue(unsigned int elem_id) const;
  virtual Real getMappedElementalValue(unsigned int elem_id) const;

protected:
  void readFile();

  const MaterialProperty<RankTwoTensor> & _mat_prop;
  const unsigned int _index_i;
  const unsigned int _index_j;

  std::vector<Real> _elem_integrals;

  FileName _file_name;
  std::map<unsigned int, std::vector<unsigned int>> _mapped_elem_id;
};

#endif // GOLEMRANKTWOMATERIALPROPERTY_H
