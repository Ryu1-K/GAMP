function y = MMSEdet(RX,CH)
y = CH.H'*inv(CH.H*CH.H'+CH.N0*eye(SIM.N))*RX.Y;
