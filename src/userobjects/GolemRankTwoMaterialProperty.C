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
#include "GolemRankTwoMaterialProperty.h"
#include "MooseMesh.h"

registerMooseObject("GolemApp", GolemRankTwoMaterialProperty);

template <>
InputParameters
validParams<GolemRankTwoMaterialProperty>()
{
  InputParameters params = validParams<ElementIntegralUserObject>();
  params.addClassDescription("UserObject that provide the value of a property from a RankTwoTensor "
                             "and mapped onto an element in case by reading a map create by "
                             "GolemMap User Object and stored in the file_name file.");
  params.addRequiredParam<MaterialPropertyName>(
      "mat_prop", "The name of the RankTWoTensor material properties to be used.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i",
      "index_i >= 0 & index_i <= 2",
      "The index i of ij for the RankTwoTensor material (0, 1, 2)");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j",
      "index_j >= 0 & index_j <= 2",
      "The index j of ij for the RankTwoTensor material (0, 1, 2)");
  params.addRequiredParam<FileName>("file_name",
                                    "The name of the file which stores the relative mapping");
  return params;
}

GolemRankTwoMaterialProperty::GolemRankTwoMaterialProperty(const InputParameters & parameters)
  : ElementIntegralUserObject(parameters),
    _mat_prop(getMaterialProperty<RankTwoTensor>("mat_prop")),
    _index_i(getParam<unsigned int>("index_i")),
    _index_j(getParam<unsigned int>("index_j")),
    _file_name(getParam<FileName>("file_name"))
{
}

void
GolemRankTwoMaterialProperty::readFile()
{
  std::string line;
  MooseUtils::checkFileReadable(_file_name);
  std::ifstream stream(_file_name);
  if (!stream.good())
    mooseError(name(), ": error opening file '", _file_name);
  while (std::getline(stream, line))
  {
    std::istringstream iss(line);
    std::vector<unsigned int> tmp;
    unsigned int ii;
    while (iss >> ii)
      tmp.push_back(ii);
    ii = tmp.front();
    tmp.erase(tmp.begin());
    _mapped_elem_id.insert(std::make_pair(ii, tmp));
    tmp.clear();
  }
  stream.close();
}

void
GolemRankTwoMaterialProperty::initialize()
{
  ElementIntegralUserObject::initialize();
  _elem_integrals.clear();
  _elem_integrals.resize(_subproblem.mesh().getMesh().max_elem_id());
  _mapped_elem_id.clear();
  readFile();
}

void
GolemRankTwoMaterialProperty::execute()
{
  Real integral_value = computeQpIntegral();
  _elem_integrals[_current_elem->id()] = integral_value;
}

void
GolemRankTwoMaterialProperty::finalize()
{
  gatherSum(_elem_integrals);
}

void
GolemRankTwoMaterialProperty::threadJoin(const UserObject & y)
{
  ElementIntegralUserObject::threadJoin(y);
  const GolemRankTwoMaterialProperty & mat_uo =
      dynamic_cast<const GolemRankTwoMaterialProperty &>(y);
  for (unsigned int i = 0; i < _elem_integrals.size(); i++)
    _elem_integrals[i] += mat_uo._elem_integrals[i];
}

Real
GolemRankTwoMaterialProperty::computeIntegral()
{
  mooseError(name(), ": you should not enter here!");
}

Real
GolemRankTwoMaterialProperty::computeQpIntegral()
{
  return _mat_prop[_qp](_index_i, _index_j);
}

Real
GolemRankTwoMaterialProperty::getElementalValue(unsigned int elem_id) const
{
  if (_elem_integrals.size() > 0)
    return _elem_integrals[elem_id];
  else
    return 0.0;
}

Real
GolemRankTwoMaterialProperty::getMappedElementalValue(unsigned int elem_id) const
{
  if (_elem_integrals.size() <= 0)
    return 0.0;
  Real value = 0.0;
  if (_mapped_elem_id.find(elem_id) == _mapped_elem_id.end())
    mooseError(name(), " Could not find the key ", elem_id, " in the current map");
  std::vector<unsigned int> map = _mapped_elem_id.find(elem_id)->second;
  for (unsigned i = 0; i < map.size(); i++)
    value += getElementalValue(map[i]);
  return value / map.size();
}
