#!/bin/sh

function compose() {
	echo "static const char *fragment_shaders[Prog_COUNT+2] = {"
	for f; do
		cat $f | sed \
			-e 's/\/\/.*//' \
			-e 's/\\/\\\\/g' \
			-e 's/^\s*//' \
			-e 's/\s\{2,\}//g' \
			-e 's/\s*\([{}();*+=,.-]\+\)\s*/\1/g' \
			-e 's/^/"/' \
			-e 's/\([a-zA-Z0-9\\]\)$/\1\\n/' \
			-e 's/$/"/'
			#-e ':a;N;$!ba;s/\n/ /g' \

		echo ", "
	done
	echo "};"
}

compose generator.glsl eroder.glsl planner.glsl tracer.glsl postprocessor.glsl common.glsl > shaders.h
