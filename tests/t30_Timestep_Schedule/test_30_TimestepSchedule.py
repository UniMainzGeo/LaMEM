import os

import pyTestHarness.unittest as pth

import pyTestHarness.launcher as launch

import re


def test_TS_Schedule():
  ranks = 4

  launch = '../bin/opt/LaMEM -ParamFile ./t30_Timestep_Schedule/TS_Schedule.dat' # This must be a relative path with respect to runLaMEM_Tests.p

  expected_file = 't30_Timestep_Schedule/t30_TS_Schedule.expected'


  def comparefunc(unittest):

    key = re.escape("Actual time step :")

    unittest.compareFloatingPoint(key,1e6)


  # Create unit test object

  ex1 = pth.pthUnitTest('t30_schedule',ranks,launch,expected_file)

  ex1.setVerifyMethod(comparefunc)

  ex1.appendKeywords('@')


  return(ex1)