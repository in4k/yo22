#!/bin/sh

function program() {
	echo "static const unsigned char prg[] = {"
	echo "CMD(0, TexNoise8U, SamplerBinding_Noise),"
  	echo "CMD(0, TexTerrain0, SamplerBinding_Terrain),"
 	echo "CMD(0, TexPlan, SamplerBinding_Plan),"
	echo "CMD(0, TexFrame, SamplerBinding_Frame),"
	echo "CMD(1, TexTerrain0, 0),"
    echo "CMD(2, ProgTerrainGenerate, 0),"
    echo "CMD(3, 0, 0),"

	#for i in `gseq 0 1024`
	for ((i = 0 ; i < 1024 ; i++ ))
	do
	  echo "CMD(1, TexTerrain1, 0),"
      echo "CMD(2, ProgTerrainErode, 0),"
      echo "CMD(0, TexTerrain0, SamplerBinding_Terrain),"
      echo "CMD(3, 0, 0),"
      echo "CMD(1, TexTerrain0, 0),"
      echo "CMD(2, ProgTerrainErode, 0),"
      echo "CMD(0, TexTerrain1, SamplerBinding_Terrain),"
      echo "CMD(3, 0, 0),"
	done

	echo "};"
}

function compose() {
	echo "static const char *fragment_shaders[Prog_COUNT+1] = {"
	for f; do
		cat $f | sed \
			-e 's/\/\/.*//' \
			-e 's/\\/\\\\/g' \
			-e 's/^\s*//' \
			-e 's/\s\{2,\}//g' \
			-e 's/\s*\([{}(;)*+=,.-]\+\)\s*/\1/g' \
			-e 's/^/"/' \
			-e 's/\([a-zA-Z0-9\\]\)$/\1 /' \
			-e 's/#\(.*\)/\\n#\1\\n/' \
			-e 's/$/"/'

			#-e 's/float/F/g' -e 's/vec2/V2/g' -e 's/vec3/V3/g' -e 's/vec4/V4/g' \
			#-e 's/f_loat/float/g' -e 's/v_ec2/vec2/g' -e 's/v_ec3/vec3/g' -e 's/v_ec4/vec4/g' \

			#-e ':a;N;$!ba;s/\n/ /g' \

		echo ", "
	done
	echo "};"
}

compose generator.glsl NOeroder.glsl planner.glsl tracer.glsl postprocessor.glsl common.glsl > shaders.h
