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
#include "GolemMatProp.h"

registerMooseObject("GolemApp", GolemMatProp);

template <>
InputParameters
validParams<GolemMatProp>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription(
      "Aux Kernel that stores the value of a material also for mapped elements.");
  params.addRequiredParam<UserObjectName>("material_user_object",
                                          "The user object to retrieve the material from.");
  return params;
}

GolemMatProp::GolemMatProp(const InputParameters & parameters)
  : AuxKernel(parameters),
    _mat_uo(getUserObject<GolemRankTwoMaterialProperty>("material_user_object"))
{
}

Real
GolemMatProp::computeValue()
{
  Real value = 0.0;
  if (_current_elem->dim() == _mesh.dimension())
    value = _mat_uo.getElementalValue(_current_elem->id());
  else if (_current_elem->dim() == (_mesh.dimension() - 1))
    value = _mat_uo.getMappedElementalValue(_current_elem->id());

  return value;
}
