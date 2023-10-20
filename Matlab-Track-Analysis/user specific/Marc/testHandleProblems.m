function testHandleProblems(eset)
tic;
%{
hs = eset.gatherField('headSwing');
toc
mean([hs.sign])
%}
mean (eset.gatherSubField('headSwing', 'sign'))
toc
