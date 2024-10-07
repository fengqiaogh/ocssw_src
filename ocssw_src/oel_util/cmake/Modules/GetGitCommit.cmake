# Get Git commit hash for head branch
#
# get_git_commit(<result> [SHORT [<length>]])

function(get_git_commit result)
	set(out "")
	if (ARGV1 STREQUAL SHORT)
		if (ARGV2)
			execute_process(COMMAND git rev-parse --short=${ARGV2} HEAD OUTPUT_VARIABLE out)
		else ()
			execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE out)
		endif ()
	else ()
		execute_process(COMMAND git rev-parse HEAD OUTPUT_VARIABLE out)
	endif ()
	string(STRIP ${out} out)
	set(${result} ${out} PARENT_SCOPE)
endfunction(get_git_commit)

