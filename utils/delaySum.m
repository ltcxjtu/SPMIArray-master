function Sf = delaySum(ftbin,Df)

Sf = squeeze(sum(bsxfun(@times,ftbin,Df),1));

end
