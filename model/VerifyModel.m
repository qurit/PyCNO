total = 0;
LambdaPhys = 0.00632;
nmol2MBQ = LambdaPhys / 60 * 6.022e23 / 10^9 / 10^6;
figure
hold on
count = 0;
for i = 1:length(results.DataInfo)-20
if contains(results.DataInfo{i, 1}.Name, 'Hot')
total = total + results.Data(:,i);
plot(results.Time, results.Data(:,i)*nmol2MBQ);
count = count + 1;
end
end
disp(count)
plot(results.Time, total*nmol2MBQ)
plot(results.Time, 318.4*exp(-results.Time*LambdaPhys))