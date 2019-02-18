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
#include "GolemMap.h"
#include "MooseMesh.h"

#include <fstream>

registerMooseObject("GolemApp", GolemMap);

template <>
InputParameters
validParams<GolemMap>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addClassDescription(
      "User Object that creates a map for the common ids between paired elements and write the map "
      "into a csv file to be used by GolemRankTwoMaterialProperty User Object.");
  params.addRequiredParam<FileName>("file_name", "The name of the file where to write the map.");
  params.addParam<Real>("tolerance", 1e-2, "The tolerance for the searching routine");
  return params;
}

GolemMap::GolemMap(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _file_name(getParam<FileName>("file_name")),
    _tolerance(getParam<Real>("tolerance"))
{
}

void
GolemMap::initialSetup()
{
  _matrix_elem_id.clear();
  _frac_elem_id.clear();
  for (const auto & elem : _subproblem.mesh().getMesh().element_ptr_range())
  {
    bool found_elem = (elem != nullptr);
    this->comm().max(found_elem);
    if (!found_elem)
      mooseError(name(), ": the given element ", elem->id(), " could not be found.");
    if (elem->dim() == _subproblem.mesh().getMesh().mesh_dimension())
      _matrix_elem_id.push_back(elem->id());
    else if (elem->dim() == (_subproblem.mesh().getMesh().mesh_dimension() - 1))
      _frac_elem_id.push_back(elem->id());
  }
  _mapped_elem_id.clear();
  for (auto & fel_id : _frac_elem_id)
  {
    if (!mapElements(fel_id))
      mooseError(name(), ": the element ", fel_id, " does not share any common element");
  }
  writeFile();
  _mapped_elem_id.clear();
}

bool
GolemMap::mapElements(const unsigned int id)
{
  auto fel = _subproblem.mesh().getMesh().query_elem_ptr(id);
  std::vector<unsigned int> fill_list;
  for (auto mel_id : _matrix_elem_id)
  {
    unsigned int counter = 0;
    auto mel = _subproblem.mesh().getMesh().query_elem_ptr(mel_id);
    for (unsigned int i = 0; i < fel->n_nodes(); i++)
      for (unsigned int j = 0; j < mel->n_nodes(); j++)
        if (pointsAreEqual(fel->point(i), mel->point(j), _tolerance))
          counter++;
    if (counter == fel->n_nodes())
      fill_list.push_back(mel_id);
  }
  if (fill_list.size() > 0)
  {
    _mapped_elem_id.insert(std::make_pair(id, fill_list));
    fill_list.clear();
    return true;
  }

  return false;
}

bool
GolemMap::pointsAreEqual(const Point pt0, const Point pt1, const Real tol)
{
  return (MooseUtils::absoluteFuzzyEqual(pt0(0), pt1(0), tol) &&
          MooseUtils::absoluteFuzzyEqual(pt0(1), pt1(1), tol) &&
          MooseUtils::absoluteFuzzyEqual(pt0(2), pt1(2), tol));
}

void
GolemMap::writeFile()
{
  FILE * file = fopen(_file_name.c_str(), "w");
  for (unsigned i = 0; i < _frac_elem_id.size(); i++)
  {
    std::vector<std::string> ids;
    std::stringstream ss;
    ss << _mapped_elem_id.find(_frac_elem_id[i])->first << " ";
    std::vector<unsigned int> second = _mapped_elem_id.find(_frac_elem_id[i])->second;
    for (unsigned int i = 0; i < second.size(); i++)
      ss << second[i] << " ";
    ids.push_back(ss.str());
    for (unsigned int i = 0; i < ids.size(); i++)
      fputs(ids[i].c_str(), file);
    fputs("\n", file);
  }
  fclose(file);
}
