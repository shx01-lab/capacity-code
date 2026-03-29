dA = 2; dB = 2;
theta = 0.01;
theta_2=0.004;
gamma = (0.5:0.01:0.99);
% 示例：随机Hermitian M（你用已知的AB系统矩阵替换即可）
n = dA*dB;
m = length(gamma);

% 方法A：使用3D数组存储所有矩阵（专业做法）
M = zeros(4, 4, m);
N = zeros(4, 4, m);% 创建4×4×9801的3D数组
ts = zeros(1, m);
rs = zeros(1,m);
statuses = cell(1, m);
for idx = 1:m
    M(:, :, idx) = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]+gamma(idx)* [0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 -1]+sqrt(1-gamma(idx))*[0 0 0 1;0 0 0 0;0 0 0 0; 1 0 0 0];
     try
        [ts(idx)] = sollve_AB_SDP(M(:, :, idx), theta, dA, dB);
    catch ME
        fprintf('gamma=%.4f 时出错: %s\n', gamma(idx), ME.message);
        ts(idx) = NaN;  % 标记为无效值
        statuses{idx} = ['ERROR: ', ME.message];
    end
end
for idx = 1:m
    N(:, :, idx) = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]+gamma(idx)* [0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 -1]+sqrt(1-gamma(idx))*[0 0 0 1;0 0 0 0;0 0 0 0; 1 0 0 0];
     try
        [rs(idx)] = sollve_AB_SDP(N(:, :, idx), theta_2, dA, dB);
    catch ME
        fprintf('gamma=%.4f 时出错: %s\n', gamma(idx), ME.message);
        rs(idx) = NaN;  % 标记为无效值
        statuses{idx} = ['ERROR: ', ME.message];
    end
end
plot(gamma, log2(ts), 'r', gamma, log2(rs), 'c');
legend('theta=0.01','theta=0.004')
xlabel('$\gamma$', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('c(N)', 'FontSize', 12);
cvx_clear
