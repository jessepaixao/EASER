function  fobj = fobjHinf(x,model)

model.param.set('Hb', [num2str(x,10),'[mm]']');
model.study('std3').run;
fn=mphglobal(model,'abs(freq)','dataset','dset5');

model.study('std3').feature('frmod').set('plist', fn(1:2));
model.study('std3').run;

data=mphplot(model,'pg7','createplot','off');

h=data{1, 1}{1, 1}.d;
    
fobj=double(norm(h,inf));


