/**
 * \file state.h
 */


#ifndef STATE_H_
#define STATE_H_

namespace Optimization
{
  namespace Local
  {
    /**
     * Optimization state for an optimizer
     */
    template<class DataType_, class ParameterType_>
    struct State
    {
      typedef DataType_ DataType;
      typedef ParameterType_ ParameterType;

      ParameterType best_parameters;
      DataType best_value;

      ParameterType current_parameters;
      DataType current_value;

      ParameterType former_parameters;
      DataType former_value;

      long iteration;
    };
  }
}


#endif /* STATE_H_ */
