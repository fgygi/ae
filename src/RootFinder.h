////////////////////////////////////////////////////////////////////////////////
//
// RootFinder.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: RootFinder.h,v 1.2 2009-02-14 22:47:30 fgygi Exp $

#ifndef ROOTFINDER_H
#define ROOTFINDER_H

class RootFinder
{
  private:

  bool hunting_phase, know_previous_f;
  double step, previous_x, previous_f;

  public:

  RootFinder(double initial_step)
  {
    reset(initial_step);
  }

  void reset(double initial_step)
  {
    hunting_phase = true;
    know_previous_f = false;
    step = initial_step;
  }


  double next(double x, double f)
  {
    if ( know_previous_f )
    {
      if ( f * previous_f < 0.0 )
      {
        // root found between x and previous_x
        // switch to bisection phase
        hunting_phase = false;
        step *= -0.5;
        // make step in the direction of the smaller absolute value of f
        previous_x = x;
        previous_f = f;
        return x + step;
      }
      else
      {
        if ( hunting_phase )
        {
          step *= 1.2;
          if ( f*f > previous_f*previous_f )
            step = - step;
        }
        previous_x = x;
        previous_f = f;
        return x + step;
      }
    }
    else
    {
      know_previous_f = true;
      previous_x = x;
      previous_f = f;
      return x + step;
    }
  }
};
#endif
