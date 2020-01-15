% plot maxrandomx for different lengths

len = 15;

maxRand = nan(1,len);

t = 2:len+1;

for k = t
  maxRand(k-1) = findMaxRandomX('ABC', k, 3, 0.75, 0.35);
end

figure;

plot(t, maxRand, 'r-', 'LineWidth', 2.5)
ylabel('max random x')
axis equal