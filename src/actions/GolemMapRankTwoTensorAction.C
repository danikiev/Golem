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
#include "GolemMapRankTwoTensorAction.h"
#include "Factory.h"
#include "FEProblem.h"
#include "MooseMesh.h"

#include "libmesh/string_to_enum.h"

registerMooseAction("GolemApp", GolemMapRankTwoTensorAction, "add_aux_variable");
registerMooseAction("GolemApp", GolemMapRankTwoTensorAction, "add_aux_kernel");
registerMooseAction("GolemApp", GolemMapRankTwoTensorAction, "add_user_object");

template <>
InputParameters
validParams<GolemMapRankTwoTensorAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::vector<VariableName>>(
      "rank_two_material_property",
      "The name of the (RankTwoTensor) material property which you want to map.");
  params.addRequiredParam<std::vector<unsigned int>>(
      "index_i", "The list for the first index of the RankTwoTensor.");
  params.addRequiredParam<std::vector<unsigned int>>(
      "index_j", "The list for the second index of the RankTwoTensor.");
  params.addRequiredParam<FileName>("file_name",
                                    "The name of the file where to store the mapping.");
  params.addParam<bool>("create_map",
                        true,
                        "Whether or not to create the map. If false, it assumes it is provided in "
                        "the file file_name.");
  return params;
}

GolemMapRankTwoTensorAction::GolemMapRankTwoTensorAction(InputParameters params)
  : Action(params),
    _rank_two_material_property(getParam<std::vector<VariableName>>("rank_two_material_property")),
    _index_i(getParam<std::vector<unsigned int>>("index_i")),
    _index_j(getParam<std::vector<unsigned int>>("index_j")),
    _file_name(getParam<FileName>("file_name")),
    _create_map(getParam<bool>("create_map"))
{
  if (_index_i.size() != _index_j.size())
    mooseError(name(), ": the size of the list of indices (index_i and index_j) is not the same.");
  std::string sx = "x";
  std::string sy = "y";
  std::string sz = "z";
  for (unsigned i = 0; i < _index_i.size(); i++)
  {
    if (_index_i[i] == 0)
      _reindex_i.push_back(sx);
    else if (_index_i[i] == 1)
      _reindex_i.push_back(sy);
    else if (_index_i[i] == 2)
      _reindex_i.push_back(sz);
    if (_index_j[i] == 0)
      _reindex_j.push_back(sx);
    else if (_index_j[i] == 1)
      _reindex_j.push_back(sy);
    else if (_index_j[i] == 2)
      _reindex_j.push_back(sz);
  }
}

void
GolemMapRankTwoTensorAction::act()
{
  if (_current_task == "add_aux_variable")
    createAuxVariableActions();
  else if (_current_task == "add_aux_kernel")
    createAuxKernelActions();
  else if (_current_task == "add_user_object")
    createUserObjectActions();
}

void
GolemMapRankTwoTensorAction::createAuxVariableActions()
{
  for (unsigned i = 0; i < _rank_two_material_property.size(); i++)
    for (unsigned j = 0; j < _index_i.size(); j++)
    {
      std::string name = _rank_two_material_property[i] + "_" + _reindex_i[j] + _reindex_j[j];
      _aux_variables.push_back(name);
      _problem->addAuxVariable(_aux_variables.back(),
                               FEType(Utility::string_to_enum<Order>("CONSTANT"),
                                      Utility::string_to_enum<FEFamily>("MONOMIAL")));
    }
}

void
GolemMapRankTwoTensorAction::createAuxKernelActions()
{
  std::string type = "GolemMatProp";
  InputParameters params = emptyInputParameters();
  for (unsigned i = 0; i < _rank_two_material_property.size(); i++)
    for (unsigned j = 0; j < _index_i.size(); j++)
    {
      std::string name = _rank_two_material_property[i] + "_" + _reindex_i[j] + _reindex_j[j];
      _aux_kernels.push_back(name + "_aux");
      params = _factory.getValidParams(type);
      params.set<AuxVariableName>("variable") = name;
      params.set<UserObjectName>("material_user_object") = name + "_uo";
      params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
      _problem->addAuxKernel(type, _aux_kernels.back(), params);
    }
}

void
GolemMapRankTwoTensorAction::createUserObjectActions()
{
  InputParameters params = emptyInputParameters();
  std::string type;
  if (_create_map)
  {
    type = "GolemMap";
    _user_objects.push_back("map_uo");
    params = _factory.getValidParams(type);
    params.set<FileName>("file_name") = _file_name;
    params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;
    _problem->addUserObject(type, _user_objects.back(), params);
    type = "GolemRankTwoMaterialProperty";
  }

  type = "GolemRankTwoMaterialProperty";
  for (unsigned i = 0; i < _rank_two_material_property.size(); i++)
    for (unsigned j = 0; j < _index_i.size(); j++)
    {
      std::string name = _rank_two_material_property[i] + "_" + _reindex_i[j] + _reindex_j[j];
      _user_objects.push_back(name + "_uo");
      params = _factory.getValidParams(type);
      params.set<MaterialPropertyName>("mat_prop") = _rank_two_material_property[i];
      params.set<unsigned int>("index_i") = _index_i[j];
      params.set<unsigned int>("index_j") = _index_j[j];
      params.set<FileName>("file_name") = _file_name;
      params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
      _problem->addUserObject(type, _user_objects.back(), params);
    }
}
