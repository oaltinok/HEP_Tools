#==========================
# run500
#==========================
find /minerva/data/users/oaltinok/PIDStudies/Outliers_NotRemoved/grid/central_value/minerva/ana/v10r6p12/00/00/05/00 -name "*Ana_PID*.root" -printf "%p\n" | sort > pl_run500_Outliers_NotRemoved.dat
find /minerva/data/users/oaltinok/PIDStudies/Outliers_Tammy/grid/central_value/minerva/ana/v10r6p12/00/00/05/00 -name "*Ana_PID*.root" -printf "%p\n" | sort > pl_run500_Outliers_Tammy.dat
find /minerva/data/users/oaltinok/PIDStudies/Outliers_Brandon/grid/central_value/minerva/ana/v10r6p12/00/00/05/00 -name "*Ana_PID*.root" -printf "%p\n" | sort > pl_run500_Outliers_Brandon.dat
find /minerva/data/users/oaltinok/PIDStudies/Outliers_Ozgur/grid/central_value/minerva/ana/v10r6p12/00/00/05/00 -name "*Ana_PID*.root" -printf "%p\n" | sort > pl_run500_Outliers_Ozgur.dat

#==========================
# run501
#==========================
find /minerva/data/users/oaltinok/PIDStudies/Outliers_NotRemoved/grid/central_value/minerva/ana/v10r6p12/00/00/05/01 -name "*Ana_PID*.root" -printf "%p\n" | sort > pl_run501_Outliers_NotRemoved.dat
find /minerva/data/users/oaltinok/PIDStudies/Outliers_Tammy/grid/central_value/minerva/ana/v10r6p12/00/00/05/01 -name "*Ana_PID*.root" -printf "%p\n" | sort > pl_run501_Outliers_Tammy.dat
find /minerva/data/users/oaltinok/PIDStudies/Outliers_Brandon/grid/central_value/minerva/ana/v10r6p12/00/00/05/00 -name "*Ana_PID*.root" -printf "%p\n" | sort > pl_run501_Outliers_Brandon.dat
find /minerva/data/users/oaltinok/PIDStudies/Outliers_Ozgur/grid/central_value/minerva/ana/v10r6p12/00/00/05/01 -name "*Ana_PID*.root" -printf "%p\n" | sort > pl_run501_Outliers_Ozgur.dat
