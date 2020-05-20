nSessions = numel(GLMData.Sessions);
areas = {};
for s=1:nSessions % get rid of areas without valve onset data
	for a=1:numel(GLMData.Sessions(s).Areas)
		area = GLMData.Sessions(s).Areas(a).area;
		if sum(strcmp(areas, GLMData.Sessions(s).Areas(a).area))==0
			areas = cat(1,areas, GLMData.Sessions(s).Areas(a).area);
			sIndices.(area) = []; 
			aIndices.(area) = [];
		end
		sIndices.(area) = cat(1, sIndices.(area), s);
		aIndices.(area) = cat(1, aIndices.(area), a);
	end
end

for i=1:numel(areas)
	area = areas(i);
	disp(area)
	GLMDataArea.Sessions = GLMData.Sessions(sIndices.(area{1}));
	for s=1:numel(GLMDataArea.Sessions)
		GLMDataArea.Sessions(s).Areas = GLMDataArea.Sessions(s).Areas(aIndices.(area{1})(s));
	end
	save(strcat('/lcncluster/muscinel/GDplusMemory/GLMDataV9',area{1}),'GLMDataArea','-v7.3');
end
