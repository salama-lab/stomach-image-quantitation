function thresvalue=minminT(I)

thresvalue = max([min(max(I,[],1)) min(max(I,[],2))])	 ;