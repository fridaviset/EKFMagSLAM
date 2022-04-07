figure; clf;
for a=0.3:0.1:0.8
plot(1:10,a.*(1:10),'Color',a.*ones(3,1));
hold on;
end