function index = AlterSign(x)

for i = 1:length(x)
    if x(i) > 0 && x(i+1) <= 0
        index = i+1;
        break;
    end
end