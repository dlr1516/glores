/**
 * \file criteria.h
 */

#ifndef CRITERIA_H_
#define CRITERIA_H_

#include <cmath>

namespace Optimization
{
  namespace Local
  {
    struct IterationCriterion
    {
    private:
      int max_iterations;
    public:
      IterationCriterion(long max_iterations)
      :max_iterations(max_iterations)
      {
      }

      template<class State>
      bool operator()(const State& state)
      {
        bool result = state.iteration < max_iterations;
#if VERBOSE > 5
        std::cout << "IterationCriteria: " << result << std::endl;
#endif
        return result;
      }
    };

    template<class DataType>
    struct RelativeValueCriterion
    {
    private:
      DataType ftol;
    public:
      RelativeValueCriterion(DataType ftol) :
          ftol(ftol)
      {
      }

      template<class State>
      bool operator()(const State& state)
      {
        bool result = 2 * std::abs(state.current_value - state.former_value)
        / (std::abs(state.current_value) + std::abs(state.former_value) + std::numeric_limits<DataType>::epsilon()) > ftol;
#if VERBOSE > 5
        std::cout << "RelativeValueCriterion: " << result << std::endl;
#endif
        return result;
      }
    };

    template<class Criteria1, class Criteria2>
    struct AndCriteria
    {
      Criteria1 criteria1;
      Criteria2 criteria2;

      AndCriteria(const Criteria1& criteria1, const Criteria2& criteria2) :
          criteria1(criteria1), criteria2(criteria2)
      {
      }

      template<class State>
      bool operator()(const State& state)
      {
        return criteria1(state) && criteria2(state);
      }
    };

    template<class Criteria1, class Criteria2>
    AndCriteria<Criteria1, Criteria2> make_and_criteria(const Criteria1& criteria1, const Criteria2& criteria2)
    {
      return AndCriteria<Criteria1, Criteria2>(criteria1, criteria2);
    }
  }
}

#endif /* CRITERIA_H_ */
