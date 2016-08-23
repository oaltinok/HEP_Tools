function [ ] = WeightsFromSideBands(  )

x = [1,2,3,4,5,6,7,8];
SinglePiPlus = [1,1.2,1.27,1.18,1.1,1.24,1.02,1.01];
QELike = [0.99,0.63,0.78,0.57,0.63,0.5,0.61,1.03];
WithPi0 = [1,0.94,0.89,0.94,0.97,0.9,1.06,1.02];

err_SinglePiPlus = [0.05,0.1,0.09,0.1,0.08,0.06,0.08,0.09];
err_QELike = [0.1,0.15,0.15,0.15,0.11,0.07,0.11,0.15];
err_WithPi0 = [0.03,0.05,0.05,0.06,0.05,0.05,0.04,0.05];



errorbar(SinglePiPlus,err_SinglePiPlus, 'ko');
hold on;
errorbar(QELike,err_QELike, 'ro');
errorbar(WithPi0,err_WithPi0, 'bo');



end

