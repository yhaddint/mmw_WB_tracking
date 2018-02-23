clear;clc;
N = 8;
K = 2;
H = randn(N) + 1j*randn(N);
Tu = dftmtx(N);
for ii=1:5
    V = H*Tu(:,1:K);
    [Q,R] = qr(V);
    Td = Q(:,1:K);
    U = H'*Td;
    [Q,R] = qr(U);
    error(ii) = norm(Q(:,1:K) - Tu(:,1:K),'fro')/norm(Tu(:,1:K),'fro');
    Tu = Q(:,1:K);
end
[U,S,V] = svd(H);
Td(:,1:K)'*H*Tu(:,1:K)

U(:,1:K)'*H*V(:,1:K)


figure
plot(error);
grid on
xlabel('Iteration Number')
ylabel('Difference in Updating')