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
#ifndef GOLEMMAPRANKTWOTENSORACTION_H
#define GOLEMMAPRANKTWOTENSORACTION_H

#include "Action.h"

class GolemMapRankTwoTensorAction;

template <>
InputParameters validParams<GolemMapRankTwoTensorAction>();

class GolemMapRankTwoTensorAction : public Action
{
public:
  GolemMapRankTwoTensorAction(InputParameters params);
  virtual void act() override;

protected:
  virtual void createAuxVariableActions();
  virtual void createAuxKernelActions();
  virtual void createUserObjectActions();

  std::vector<VariableName> _rank_two_material_property;
  std::vector<unsigned int> _index_i;
  std::vector<unsigned int> _index_j;
  std::vector<std::string> _reindex_i;
  std::vector<std::string> _reindex_j;

  FileName _file_name;
  bool _create_map;

  std::vector<std::string> _aux_variables;
  std::vector<std::string> _aux_kernels;
  std::vector<std::string> _user_objects;
};

#endif // GOLEMMAPRANKTWOTENSORACTION_H
