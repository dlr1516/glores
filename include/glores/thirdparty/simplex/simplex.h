/**
 * \file simplex.h
 */

#ifndef SIMPLEX_H_
#define SIMPLEX_H_

#include "state.h"

namespace Optimization
{
  namespace Local
  {
    /**
     * The Nelder-Mead Simplex algorithm
     * Complies to a simple and standard interface
     *
     * @param DataType is the type of inner values to consider
     * @param ParameterType is the type of ParameterType (Eigen if possible)
     * @param Function is the type of Function to optimize
     * @param Criterion is the type of the stopping Criterion
     */
    template<class DataType, class ParameterType, class Function, class Criterion>
    class Simplex
    {
      State<DataType, ParameterType> state;
      DataType delta;
      ParameterType deltas;
      bool use_deltas;
      Eigen::Array<DataType, Eigen::Dynamic, Eigen::Dynamic> polytope_points;
      Eigen::Array<DataType, Eigen::Dynamic, Eigen::Dynamic> polytope_values;

      Criterion criterion;

      void initialize_polytope(const ParameterType& start_point, DataType delta, const Function& fun)
      {
        polytope_points.resize(start_point.size(), start_point.size() + 1);
        polytope_values.resize(1, start_point.size() + 1);
        polytope_points.col(0) = start_point;
        polytope_values(0, 0) = fun(start_point);
        display(fun, polytope_points.col(0));
        for (int i = 1; i < start_point.size() + 1; ++i)
        {
          polytope_points.col(i) = start_point;
          polytope_points(i - 1, i) += delta;
          polytope_values(0, i) = fun(polytope_points.col(i));
          display(fun, polytope_points.col(i));
        }
      }

      void initialize_polytope(const ParameterType& start_point, ParameterType deltas, const Function& fun)
      {
        polytope_points.resize(start_point.size(), start_point.size() + 1);
        polytope_values.resize(1, start_point.size() + 1);
        polytope_points.col(0) = start_point;
        polytope_values(0, 0) = fun(start_point);
        display(fun, polytope_points.col(0));
        for (int i = 1; i < start_point.size() + 1; ++i)
        {
          polytope_points.col(i) = start_point;
          polytope_points(i - 1, i) += deltas(i - 1);
          polytope_values(0, i) = fun(polytope_points.col(i));
          display(fun, polytope_points.col(i));
        }
      }

      void display(const Function& fun, const ParameterType& parameters)
      {
#if VERBOSE > 5
        std::cout << "Point: " << parameters << std::endl;
        std::cout << "Value: " << fun(parameters) << std::endl;
#endif
      }

      void find_best_worst_near_worst(const Function& fun, int& best, int& worst, int& near_worst)
      {
        if (polytope_values(0, 0) > polytope_values(0, 1))
        {
          worst = 0;
          near_worst = 1;
        }
        else
        {
          worst = 1;
          near_worst = 0;
        }
        best = near_worst;
        for (int i = 2; i < polytope_values.cols(); ++i)
        {
          if (polytope_values(0, i) < polytope_values(0, best))
          {
            best = i;
          }
          if (polytope_values(0, i) > polytope_values(0, worst))
          {
            near_worst = worst;
            worst = i;
          }
          else if (polytope_values(0, i) > polytope_values(0, near_worst))
          {
            near_worst = i;
          }
        }
#if VERBOSE > 5
        std::cout << "Worst value is in position " << worst << "\t" << polytope_values(0, worst) << std::endl << polytope_points.col(worst) << std::endl;
        std::cout << "Near-worst value is in position " << near_worst << "\t" << polytope_values(0, near_worst) << std::endl << polytope_points.col(near_worst) << std::endl;
        std::cout << "Best value is in position " << best << "\t" << polytope_values(0, best) << std::endl << polytope_points.col(best) << std::endl;
#endif
      }

      ParameterType create_new_parameters(const ParameterType& sum, const ParameterType& discarded_point, DataType t)
      {
        DataType fac1 = (1 - t) / sum.size();
        DataType fac2 = fac1 - t;
        return sum * fac1 - discarded_point * fac2;
      }

    public:
      Simplex(const Criterion& criterion)
      :delta(0), use_deltas(false), criterion(criterion)
      {

      }

      void optimize(const Function& fun)
      {
        state.iteration = 0;
        state.current_value = fun(state.current_parameters);
        state.former_value = std::numeric_limits<DataType>::max();
        state.best_value = std::numeric_limits<DataType>::max();

        if(use_deltas)
        {
          initialize_polytope(state.current_parameters, deltas, fun);
        }
        else
        {
          initialize_polytope(state.current_parameters, delta, fun);
        }

        while (criterion(state))
        {
#if VERBOSE > 4
          std::cout << "Starting iteration " << state.iteration << std::endl;
#endif
          int best, worst, near_worst;
          find_best_worst_near_worst(fun, best, worst, near_worst);
          state.current_value = polytope_values(0, best);
          state.current_parameters = polytope_points.col(best);
          state.former_value = polytope_values(0, worst);
          state.former_parameters = polytope_points.col(worst);

          if (state.current_value < state.best_value)
          {
            state.best_value = state.current_value;
            state.best_parameters = state.current_parameters;
          }

          ParameterType new_parameters = create_new_parameters(polytope_points.rowwise().sum(), polytope_points.col(worst), -1);
          DataType new_value = fun(new_parameters);
#if VERBOSE > 5
          std::cout << "Trying normal" << std::endl;
          display(fun, new_parameters);
#endif
          if (new_value < state.best_value)
          {
            ParameterType expansion_parameters = create_new_parameters(polytope_points.rowwise().sum(), polytope_points.col(worst), -2);
            DataType expansion_value = fun(expansion_parameters);
#if VERBOSE > 5
            std::cout << "Trying expansion" << std::endl;
            display(fun, expansion_parameters);
#endif
            if (expansion_value < state.best_value)
            {
              state.current_value = polytope_values(0, worst) = expansion_value;
              state.current_parameters = polytope_points.col(worst) = expansion_parameters;
            }
            else
            {
              state.current_value = polytope_values(0, worst) = new_value;
              state.current_parameters = polytope_points.col(worst) = new_parameters;
            }
          }
          else if (new_value > polytope_values(0, near_worst))
          {
            // New point is not better than near worst
            ParameterType contraction_parameters = create_new_parameters(polytope_points.rowwise().sum(), polytope_points.col(worst), -.5);
            DataType contraction_value = fun(contraction_parameters);
#if VERBOSE > 5
            std::cout << "Trying contraction" << std::endl;
            display(fun, contraction_parameters);
#endif
            if (contraction_value > new_value)
            {
#if VERBOSE > 5
            std::cout << "Contraction around lowest" << std::endl;
            display(fun, contraction_parameters);
#endif
              ParameterType best_parameters = polytope_points.col(best);
              //std::cout << ((polytope_points.colwise() - best_parameters.array()) / 2) << std::endl;
              polytope_points = ((polytope_points.colwise() - best_parameters.array()) / 2).colwise() + best_parameters.array();
              for(int i = 0; i < polytope_points.cols(); ++i)
              {
                polytope_values(0, i) = fun(polytope_points.col(i));
              }
            }
            else
            {
              state.current_value = polytope_values(0, worst) = contraction_value;
              state.current_parameters = polytope_points.col(worst) = contraction_parameters;
            }
          }
          else
          {
            // New point, not the best, but better than the near worst
            state.current_value = polytope_values(0, worst) = new_value;
            state.current_parameters = polytope_points.col(worst) = new_parameters;
          }

          ++state.iteration;
        }
      }

      /**
       * Retrieves the best parameters
       */
      const ParameterType& get_best_parameters() const
      {
        return state.best_parameters;
      }

      /**
       * Retrieves the best final value
       */
      const DataType get_best_value() const
      {
        return state.best_value;
      }

      void set_start_point(const ParameterType& point)
      {
        state.current_parameters = point;
      }

      void set_delta(DataType delta)
      {
        this->delta = delta;
        use_deltas = false;
      }

      void set_delta(ParameterType deltas)
      {
        this->deltas = deltas;
        use_deltas = true;
      }

      void set_iterations(long iterations)
      {
        state.max_iterations = iterations;
      }

      int get_stop_criterion() const;
    };

    template<class Function, class Criterion>
    static Simplex<typename Function::DataType, typename Function::ParameterType, Function, Criterion> build_simplex(const Function& fun, const Criterion& criterion)
    {
      return Simplex<typename Function::DataType, typename Function::ParameterType, Function, Criterion>(criterion);
    }
  }
}

#endif /* SIMPLEX_H_ */
