function [] = plotIMU(accelData, gyroData, tVec, figNumAccel, figNumGyro)
figure(figNumAccel)
subplot(3,1,1); plot(tVec, accelData(:,1))
subplot(3,1,2); plot(tVec, accelData(:,2)); 
subplot(3,1,3); plot(tVec, accelData(:,3));
sgtitle('Accelerometer data (X,Y,Z)')

figure(figNumGyro)
subplot(3,1,1); plot(tVec, gyroData(:,1))
subplot(3,1,2); plot(tVec, gyroData(:,2)); 
subplot(3,1,3); plot(tVec, gyroData(:,3)); 
sgtitle('Gyro data (X, Y, Z)')
end