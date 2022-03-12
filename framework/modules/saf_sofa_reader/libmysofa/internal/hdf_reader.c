/*
 Copyright (c) 2016, Symonics GmbH, Christian Hoene
 All rights reserved.
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:
	 (1) Redistributions of source code must retain the above copyright
	 notice, this list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright
	 notice, this list of conditions and the following disclaimer in
	 the documentation and/or other materials provided with the
	 distribution.
	 (3)The name of the author may not be used to
	 endorse or promote products derived from this software without
	 specific prior written permission.
 THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
 INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

#if defined(SAF_ENABLE_SOFA_READER_MODULE)

#ifdef _MSC_VER
# pragma warning(disable : 4101)
# pragma warning(disable : 4244)
#endif

#include <ctype.h>
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>
#include "hdf_reader.h"

/* Only external library requirement is zlib: */
#include "saf_externals.h"

/* ========================================================================== */
/*                                 HDF Reader                                 */
/* ========================================================================== */

/* checks file address.
 * NULL is an invalid address indicating a invalid field
 */
int validAddress(struct READER *reader, uint64_t address) {
	return address > 0 && address < reader->superblock.end_of_file_address;
}

/* little endian */
uint64_t readValue(struct READER *reader, int size) {
	int i, c;
	uint64_t value;
	c = fgetc(reader->fhd);
	if (c < 0)
		return 0xffffffffffffffffLL;
	value = (uint8_t)c;
	for (i = 1; i < size; i++) {
		c = fgetc(reader->fhd);
		if (c < 0)
			return 0xffffffffffffffffLL;
		value |= ((uint64_t)c) << (i * 8);
	}
	return value;
}


/* ========================================================================== */
/*                                   BTREE                                    */
/* ========================================================================== */

/*
 *
 00000370  42 54 4c 46 00 08 00 5b  01 00 00 00 2d 00 00 07  |BTLF...[....-...|
 00000380  00 00 00 f8 ea 72 15 00  c9 03 00 00 00 26 00 00  |.....r.......&..|
 00000390  14 00 00 00 32 32 7c 17  00 22 02 00 00 00 32 00  |....22|.."....2.|
 000003a0  00 0b 00 00 00 07 ef 9c  26 00 bb 01 00 00 00 46  |........&......F|
 000003b0  00 00 09 00 00 00 e5 f6  ba 26 00 45 03 00 00 00  |.........&.E....|
 000003c0  34 00 00 11 00 00 00 f6  71 f0 2e 00 a3 02 00 00  |4.......q.......|
 000003d0  00 3e 00 00 0d 00 00 00  61 36 dc 36 00 79 03 00  |.>......a6.6.y..|
 000003e0  00 00 35 00 00 12 00 00  00 97 1b 4e 45 00 88 01  |..5........NE...|
 000003f0  00 00 00 33 00 00 08 00  00 00 56 d7 d0 47 00 ae  |...3......V..G..|
 00000400  03 00 00 00 1b 00 00 13  00 00 00 2f 03 50 5a 00  |.........../.PZ.|
 00000410  22 01 00 00 00 39 00 00  06 00 00 00 b7 88 37 66  |"....9........7f|
 00000420  00 01 03 00 00 00 28 00  00 0f 00 00 00 dc aa 47  |......(........G|
 00000430  66 00 16 04 00 00 00 2c  00 00 15 00 00 00 6b 54  |f......,......kT|
 00000440  7d 77 00 fd 00 00 00 00  25 00 00 05 00 00 00 7d  |}w......%......}|
 00000450  0c 8c 9e 00 29 03 00 00  00 1c 00 00 10 00 00 00  |....)...........|
 00000460  4c f3 0e a0 00 16 00 00  00 00 25 00 00 00 00 00  |L.........%.....|
 00000470  00 e7 30 2d ab 00 01 02  00 00 00 21 00 00 0a 00  |..0-.......!....|
 00000480  00 00 35 b5 69 b0 00 e1  02 00 00 00 20 00 00 0e  |..5.i....... ...|
 00000490  00 00 00 2b c5 8b c4 00  3b 00 00 00 00 20 00 00  |...+....;.... ..|
 000004a0  01 00 00 00 09 a0 74 cc  00 93 00 00 00 00 2f 00  |......t......./.|
 000004b0  00 03 00 00 00 3f 48 ef  d6 00 5b 00 00 00 00 38  |.....?H...[....8|
 000004c0  00 00 02 00 00 00 f1 7e  7d dd 00 54 02 00 00 00  |.......~}..T....|
 000004d0  4f 00 00 0c 00 00 00 48  35 ff f5 00 c2 00 00 00  |O......H5.......|
 000004e0  00 3b 00 00 04 00 00 00  ad 61 4e ff 63 42 f7 73  |.;.......aN.cB.s|
 000004f0  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  |................|
 *
 00000570  42 54 4c 46 00 09 00 16  00 00 00 00 25 00 00 00  |BTLF........%...|
 00000580  00 00 00 00 3b 00 00 00  00 20 00 00 01 00 00 00  |....;.... ......|
 00000590  00 5b 00 00 00 00 38 00  00 02 00 00 00 00 93 00  |.[....8.........|
 000005a0  00 00 00 2f 00 00 03 00  00 00 00 c2 00 00 00 00  |.../............|
 000005b0  3b 00 00 04 00 00 00 00  fd 00 00 00 00 25 00 00  |;............%..|
 000005c0  05 00 00 00 00 22 01 00  00 00 39 00 00 06 00 00  |....."....9.....|
 000005d0  00 00 5b 01 00 00 00 2d  00 00 07 00 00 00 00 88  |..[....-........|
 000005e0  01 00 00 00 33 00 00 08  00 00 00 00 bb 01 00 00  |....3...........|
 000005f0  00 46 00 00 09 00 00 00  00 01 02 00 00 00 21 00  |.F............!.|
 00000600  00 0a 00 00 00 00 22 02  00 00 00 32 00 00 0b 00  |......"....2....|
 00000610  00 00 00 54 02 00 00 00  4f 00 00 0c 00 00 00 00  |...T....O.......|
 00000620  a3 02 00 00 00 3e 00 00  0d 00 00 00 00 e1 02 00  |.....>..........|
 00000630  00 00 20 00 00 0e 00 00  00 00 01 03 00 00 00 28  |.. ............(|
 00000640  00 00 0f 00 00 00 00 29  03 00 00 00 1c 00 00 10  |.......)........|
 00000650  00 00 00 00 45 03 00 00  00 34 00 00 11 00 00 00  |....E....4......|
 00000660  00 79 03 00 00 00 35 00  00 12 00 00 00 00 ae 03  |.y....5.........|
 00000670  00 00 00 1b 00 00 13 00  00 00 00 c9 03 00 00 00  |................|
 00000680  26 00 00 14 00 00 00 00  16 04 00 00 00 2c 00 00  |&............,..|
 00000690  15 00 00 00 d3 c7 19 a0  00 00 00 00 00 00 00 00  |................|
 000006a0  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  |................|
 */

 // LCOV_EXCL_START

static int readBTLF(struct READER *reader, struct BTREE *btree,
	int number_of_records, union RECORD *records) {

	int i;

	uint8_t type, message_flags;
	uint32_t creation_order, hash_of_name;
	uint64_t heap_id;

	char buf[5];

	UNUSED(btree);
	UNUSED(heap_id);
	UNUSED(hash_of_name);
	UNUSED(creation_order);
	UNUSED(message_flags);

	/* read signature */
	if (fread(buf, 1, 4, reader->fhd) != 4 || strncmp(buf, "BTLF", 4)) {
		mylog("cannot read signature of BTLF\n"); // LCOV_EXCL_LINE
		return MYSOFA_INVALID_FORMAT;             // LCOV_EXCL_LINE
	}
	buf[4] = 0;
	mylog("%08" PRIX64 " %.4s\n", (uint64_t)ftell(reader->fhd) - 4, buf);

	if (fgetc(reader->fhd) != 0) {
		mylog("object BTLF must have version 0\n"); // LCOV_EXCL_LINE
		return MYSOFA_INVALID_FORMAT;               // LCOV_EXCL_LINE
	}

	type = (uint8_t)fgetc(reader->fhd);

	for (i = 0; i < number_of_records; i++) {

		switch (type) {
		case 5:
			records->type5.hash_of_name = (uint32_t)readValue(reader, 4);
			records->type5.heap_id = readValue(reader, 7);
			mylog(" type5 %08X %14" PRIX64 "\n", records->type5.hash_of_name,
				records->type5.heap_id);
			records++;
			break;

		case 6:
			/*creation_order = */
			readValue(reader, 8);
			/*heap_id = */
			readValue(reader, 7);
			break;

		case 8:
			/*heap_id = */
			readValue(reader, 8);
			/*message_flags = */
			fgetc(reader->fhd);
			/*creation_order = */
			readValue(reader, 4);
			/*hash_of_name = */
			readValue(reader, 4);
			break;

		case 9:
			/*heap_id = */
			readValue(reader, 8);
			/*message_flags = */
			fgetc(reader->fhd);
			/*creation_order = */
			readValue(reader, 4);
			break;

		default:
			mylog("object BTLF has unknown type %d\n", type);
			return MYSOFA_INVALID_FORMAT;
		}
	}

	/*    fseeko(reader->fhd, bthd->root_node_address + bthd->node_size,
	 * SEEK_SET); skip checksum */

	return MYSOFA_OK;
}

/*  III.A.2. Disk Format: Level 1A2 - Version 2 B-trees
 000002d0  32 1d 42 54 48 44 00 08  00 02 00 00 11 00 00 00  |2.BTHD..........|
 000002e0  64 28 70 03 00 00 00 00  00 00 16 00 16 00 00 00  |d(p.............|
 000002f0  00 00 00 00 30 12 d9 6e  42 54 48 44 00 09 00 02  |....0..nBTHD....|
 00000300  00 00 0d 00 00 00 64 28  70 05 00 00 00 00 00 00  |......d(p.......|
 00000310  16 00 16 00 00 00 00 00  00 00 e2 0d 76 5c 46 53  |............v\FS|
 */

int btreeRead(struct READER *reader, struct BTREE *btree) {
	char buf[5];

	/* read signature */
	if (fread(buf, 1, 4, reader->fhd) != 4 || strncmp(buf, "BTHD", 4)) {
		mylog("cannot read signature of BTHD\n");
		return MYSOFA_INVALID_FORMAT;
	}
	buf[4] = 0;
	mylog("%08" PRIX64 " %.4s\n", (uint64_t)ftell(reader->fhd) - 4, buf);

	if (fgetc(reader->fhd) != 0) {
		mylog("object BTHD must have version 0\n");
		return MYSOFA_INVALID_FORMAT;
	}

	btree->type = (uint8_t)fgetc(reader->fhd);
	btree->node_size = (uint32_t)readValue(reader, 4);
	btree->record_size = (uint16_t)readValue(reader, 2);
	btree->depth = (uint16_t)readValue(reader, 2);

	btree->split_percent = (uint8_t)fgetc(reader->fhd);
	btree->merge_percent = (uint8_t)fgetc(reader->fhd);
	btree->root_node_address =
		(uint64_t)readValue(reader, reader->superblock.size_of_offsets);
	btree->number_of_records = (uint16_t)readValue(reader, 2);
	if (btree->number_of_records > 0x1000)
		return MYSOFA_UNSUPPORTED_FORMAT;
	btree->total_number =
		(uint64_t)readValue(reader, reader->superblock.size_of_lengths);

	/*    fseek(reader->fhd, 4, SEEK_CUR);  skip checksum */

	if (btree->total_number > 0x10000000)
		return MYSOFA_NO_MEMORY;
	btree->records = malloc(sizeof(btree->records[0]) * btree->total_number);
	if (!btree->records)
		return MYSOFA_NO_MEMORY;
	memset(btree->records, 0, sizeof(btree->records[0]) * btree->total_number);

	/* read records */
	if (fseek(reader->fhd, btree->root_node_address, SEEK_SET) < 0)
		return errno;
	return readBTLF(reader, btree, btree->number_of_records, btree->records);
}

// LCOV_EXCL_STOP

void btreeFree(struct BTREE *btree) { free(btree->records); }

/*  III.A.1. Disk Format: Level 1A1 - Version 1 B-trees
 *
 */

int treeRead(struct READER *reader, struct DATAOBJECT *data) {

	int i, j, err, olen, elements, size, x, y, z, b, e, dy, dz, sx, sy, sz, dzy, szy;
	char *input, *output;

	uint8_t node_type, node_level;
	uint16_t entries_used;
	uint32_t size_of_chunk;
	uint32_t filter_mask;
	uint64_t address_of_left_sibling, address_of_right_sibling,
		child_pointer, key, store;
	int start[4];
	char buf[5];

	UNUSED(node_level);
	UNUSED(address_of_right_sibling);
	UNUSED(address_of_left_sibling);
	UNUSED(key);

	if (data->ds.dimensionality > 3) {
		mylog("TREE dimensions > 3"); // LCOV_EXCL_LINE
		return MYSOFA_INVALID_FORMAT; // LCOV_EXCL_LINE
	}

	/* read signature */
	if (fread(buf, 1, 4, reader->fhd) != 4 || strncmp(buf, "TREE", 4)) {
		mylog("cannot read signature of TREE\n"); // LCOV_EXCL_LINE
		return MYSOFA_INVALID_FORMAT;             // LCOV_EXCL_LINE
	}
	buf[4] = 0;
	mylog("%08" PRIX64 " %.4s\n", (uint64_t)ftell(reader->fhd) - 4, buf);

	node_type = (uint8_t)fgetc(reader->fhd);
	node_level = (uint8_t)fgetc(reader->fhd);
	entries_used = (uint16_t)readValue(reader, 2);
	if (entries_used > 0x1000)
		return MYSOFA_UNSUPPORTED_FORMAT; // LCOV_EXCL_LINE
	address_of_left_sibling =
		readValue(reader, reader->superblock.size_of_offsets);
	address_of_right_sibling =
		readValue(reader, reader->superblock.size_of_offsets);

	elements = 1;
	for (j = 0; j < data->ds.dimensionality; j++)
		elements *= data->datalayout_chunk[j];
	dy = data->datalayout_chunk[1];
	dz = data->datalayout_chunk[2];
	sx = (int)data->ds.dimension_size[0];
	sy = (int)data->ds.dimension_size[1];
	sz = (int)data->ds.dimension_size[2];
	dzy = dz * dy;
	szy = sz * sy;
	size = data->datalayout_chunk[data->ds.dimensionality];

	mylog("elements %d size %d\n", elements, size);

	if (elements <= 0 || size <= 0 || elements >= 0x100000 || size > 0x10)
		return MYSOFA_INVALID_FORMAT; // LCOV_EXCL_LINE
	if (!(output = malloc(elements * size))) {
		return MYSOFA_NO_MEMORY; // LCOV_EXCL_LINE
	}

	for (e = 0; e < entries_used * 2; e++) {
		if (node_type == 0) {
			key = readValue(reader, reader->superblock.size_of_lengths);
		}
		else {
			size_of_chunk = (uint32_t)readValue(reader, 4);
			filter_mask = (uint32_t)readValue(reader, 4);
			if (filter_mask) {
				mylog("TREE all filters must be enabled\n"); // LCOV_EXCL_LINE
				free(output);                                // LCOV_EXCL_LINE
				return MYSOFA_INVALID_FORMAT;                // LCOV_EXCL_LINE
			}

			for (j = 0; j < data->ds.dimensionality; j++) {
				start[j] = (int)readValue(reader, 8);
				mylog("start %d %" PRIu64 "\n", j, start[j]);
			}

			if (readValue(reader, 8)) {
				break;
			}

			child_pointer = readValue(reader, reader->superblock.size_of_offsets);
			mylog(" data at %" PRIX64 " len %u\n", child_pointer, size_of_chunk);

			/* read data */
			store = ftell(reader->fhd);
			if (fseek(reader->fhd, child_pointer, SEEK_SET) < 0) {
				free(output); // LCOV_EXCL_LINE
				return errno; // LCOV_EXCL_LINE
			}

			if (!(input = malloc(size_of_chunk))) {
				free(output);            // LCOV_EXCL_LINE
				return MYSOFA_NO_MEMORY; // LCOV_EXCL_LINE
			}
			if (fread(input, 1, size_of_chunk, reader->fhd) != size_of_chunk) {
				free(output);                 // LCOV_EXCL_LINE
				free(input);                  // LCOV_EXCL_LINE
				return MYSOFA_INVALID_FORMAT; // LCOV_EXCL_LINE
			}

			olen = elements * size;
			err = gunzip(size_of_chunk, input, &olen, output);
			free(input);

			mylog("   gunzip %d %d %d\n", err, olen, elements * size);
			if (err || olen != elements * size) {
				free(output);                 // LCOV_EXCL_LINE
				return MYSOFA_INVALID_FORMAT; // LCOV_EXCL_LINE
			}

			switch (data->ds.dimensionality) {
			case 1:
				for (i = 0; i < olen; i++) {
					b = i / elements;
					x = i % elements + start[0];
					if (x < sx) {

						j = x * size + b;
						if (j >= 0 && j < data->data_len) {
							((char *)data->data)[j] = output[i];
						}
					}
				}
				break;
			case 2:
				for (i = 0; i < olen; i++) {
					b = i / elements;
					x = i % elements;
					y = x % dy + start[1];
					x = x / dy + start[0];
					if (y < sy && x < sx) {
						j = ((x * sy + y) * size) + b;
						if (j >= 0 && j < data->data_len) {
							((char *)data->data)[j] = output[i];
						}
					}
				}
				break;
			case 3:
				/* Some minor speed ups if dz==1 */
				if (dz == 1) {
					if (start[2] >= sz)
						break;
					for (i = 0; i < olen; i++) {
						b = i / elements;
						x = i % elements;
						y = x % dy + start[1];
						x = (x / dzy) + start[0];
						if (y < sy && x < sx) {
							j = (x * szy + y * sz + start[2]) * size + b;
							if (j >= 0 && j < data->data_len) {
								((char *)data->data)[j] = output[i];
							}
						}
					}
				}
				/* Some minor speed ups if dy==1 */
				else if (dy == 1) {
					if (start[1] >= sy)
						break;
					for (i = 0; i < olen; i++) {
						b = i / elements;
						x = i % elements;
						z = x % dz + start[2];
						x = (x / dzy) + start[0];
						if (z < sz  && x < sx) {
							j = (x * szy + start[1] * sz + z) * size + b;
							if (j >= 0 && j < data->data_len) {
								((char *)data->data)[j] = output[i];
							}
						}
					}
				}
				/* Otherwise, revert to the original code: */
				else {
					for (i = 0; i < olen; i++) {
						b = i / elements;
						x = i % elements;
						z = x % dz + start[2];
						y = (x / dz) % dy + start[1];
						x = (x / dzy) + start[0];
						if (z < sz && y < sy && x < sx) {
							j = (x * szy + y * sz + z) * size + b;
							if (j >= 0 && j < data->data_len) {
								((char *)data->data)[j] = output[i];
							}
						}
					}
				}
				break;
			default:
				mylog("invalid dim\n");       // LCOV_EXCL_LINE
				return MYSOFA_INTERNAL_ERROR; // LCOV_EXCL_LINE
			}

			if (fseek(reader->fhd, store, SEEK_SET) < 0) {
				free(output); // LCOV_EXCL_LINE
				return errno; // LCOV_EXCL_LINE
			}
		}
	}

	free(output);
	if (fseek(reader->fhd, 4, SEEK_CUR) < 0) /* skip checksum */
		return errno;                          // LCOV_EXCL_LINE

	return MYSOFA_OK;
}


/* ========================================================================== */
/*                                    GCOL                                    */
/* ========================================================================== */

/*  III.E. Disk Format: Level 1E - Global Heap
 */

static int readGCOL(struct READER *reader) {

	uint16_t reference_count, address;
	uint64_t collection_size, end;
	struct GCOL *gcol;
	char buf[5];

	UNUSED(reference_count);

	/* read signature */
	if (fread(buf, 1, 4, reader->fhd) != 4 || strncmp(buf, "GCOL", 4)) {
		mylog("cannot read signature of global heap collection\n");
		return MYSOFA_INVALID_FORMAT;
	}
	buf[4] = 0;

	if (fgetc(reader->fhd) != 1) {
		mylog("object GCOL must have version 1\n");
		return MYSOFA_INVALID_FORMAT;
	}
	if (fgetc(reader->fhd) < 0 || fgetc(reader->fhd) < 0 ||
		fgetc(reader->fhd) < 0)
		return MYSOFA_READ_ERROR;

	address = ftell(reader->fhd);
	end = address;
	collection_size = readValue(reader, reader->superblock.size_of_lengths);
	if (collection_size > 0x400000000) {
		mylog("collection_size is too large\n");
		return MYSOFA_INVALID_FORMAT;
	}
	end += collection_size - 8;

	while (ftell(reader->fhd) <= (long)(end - 8 - reader->superblock.size_of_lengths)) {

		gcol = malloc(sizeof(*gcol));
		if (!gcol)
			return MYSOFA_NO_MEMORY;
		gcol->heap_object_index = readValue(reader, 2);
		if (gcol->heap_object_index == 0) {
			free(gcol);
			break;
		}
		reference_count = readValue(reader, 2);
		if (fseek(reader->fhd, 4, SEEK_CUR) < 0) {
			free(gcol);
			return errno;
		}
		gcol->object_size = readValue(reader, reader->superblock.size_of_lengths);
		if (gcol->object_size > 8) {
			free(gcol);
			return MYSOFA_UNSUPPORTED_FORMAT;
		}
		gcol->value = readValue(reader, (int)gcol->object_size);
		gcol->address = address;
		mylog(" GCOL object %d size %" PRIu64 " value %08" PRIX64 "\n",
			gcol->heap_object_index, gcol->object_size, gcol->value);

		gcol->next = reader->gcol;
		reader->gcol = gcol;
	}

	mylog(" END %08lX vs. %08" PRIX64 "\n", ftell(reader->fhd),
		end); /* bug in the normal hdf5 specification */
  /*    fseek(reader->fhd, end, SEEK_SET); */
	return MYSOFA_OK;
}

int gcolRead(struct READER *reader, uint64_t gcol, int reference,
	uint64_t *dataobject) {
	long pos;
	struct GCOL *p = reader->gcol;

	while (p && p->address != gcol && p->heap_object_index != reference) {
		p = p->next;
	}
	if (!p) {
		pos = ftell(reader->fhd);
		if (fseek(reader->fhd, gcol, SEEK_SET) < 0)
			return MYSOFA_READ_ERROR;
		readGCOL(reader);
		if (pos < 0)
			return MYSOFA_READ_ERROR;
		if (fseek(reader->fhd, pos, SEEK_SET) < 0)
			return MYSOFA_READ_ERROR;

		p = reader->gcol;
		while (p && p->address != gcol && p->heap_object_index != reference) {
			p = p->next;
		}
		if (!p) {
			mylog("unknown gcol %" PRIX64 " %d\n", gcol, reference);
			return MYSOFA_INVALID_FORMAT;
		}
	}
	*dataobject = p->value;

	return MYSOFA_OK;
}
#if 0

gcol = reader->gcol;
for (;;) {
	if (gcol == NULL) {
		mylog("reference unknown!\n");
		return MYSOFA_INVALID_FORMAT;
	}
	if (gcol->heap_object_index == reference) {
		mylog("found reference at %LX\n", gcol->object_pos);
		break;
		pos = ftell(reader->fhd);
		fseek(reader->fhd, gcol->object_pos, SEEK_SET);
		dt2 = *dt;
		dt2.list = 0;
		dt2.size = gcol->object_size;
		readDataVar(reader, &dt2, ds);
		fseek(reader->fhd, pos, SEEK_SET);
		break;
	}
	gcol = gcol->next;
#endif

	void gcolFree(struct GCOL *gcol) {
		if (gcol) {
			gcolFree(gcol->next);
			free(gcol);
		}
	}


	/* ========================================================================== */
	/*                                  GUNZIP                                    */
	/* ========================================================================== */

	int gunzip(int inlen, char *in, int *outlen, char *out) {
		int err;
		z_stream stream;

		memset(&stream, 0, sizeof(stream));
		stream.avail_in = inlen;
		stream.next_in = (unsigned char *)in;
		stream.avail_out = *outlen;
		stream.next_out = (unsigned char *)out;

		err = inflateInit(&stream);
		if (err)
			return err;

        err = inflate(&stream, Z_SYNC_FLUSH);
		*outlen = (int)stream.total_out;
		inflateEnd(&stream);
		if (err && err != Z_STREAM_END) {
			mylog(" gunzip error %d %s\n", err, stream.msg);
			return err;
		}

		return MYSOFA_OK;
	}


	/* ========================================================================== */
	/*                                 SUPERBLOCK                                 */
	/* ========================================================================== */

	/* read superblock
	 00000000  89 48 44 46 0d 0a 1a 0a  02 08 08 00 00 00 00 00  |.HDF............|
	 00000010  00 00 00 00 ff ff ff ff  ff ff ff ff bc 62 12 00  |.............b..|
	 00000020  00 00 00 00 30 00 00 00  00 00 00 00 27 33 a3 16  |....0.......'3..|
	 */

	int superblockRead2or3(struct READER *reader, struct SUPERBLOCK *superblock) {
		superblock->size_of_offsets = (uint8_t)fgetc(reader->fhd);
		superblock->size_of_lengths = (uint8_t)fgetc(reader->fhd);
		if (fgetc(reader->fhd) < 0) /* File Consistency Flags */
			return MYSOFA_READ_ERROR;

		if (superblock->size_of_offsets < 2 || superblock->size_of_offsets > 8 ||
			superblock->size_of_lengths < 2 || superblock->size_of_lengths > 8) {
			mylog("size of offsets and length is invalid: %d %d\n",
				superblock->size_of_offsets, superblock->size_of_lengths);
			return MYSOFA_UNSUPPORTED_FORMAT;
		}

		superblock->base_address = readValue(reader, superblock->size_of_offsets);
		superblock->superblock_extension_address =
			readValue(reader, superblock->size_of_offsets);
		superblock->end_of_file_address =
			readValue(reader, superblock->size_of_offsets);
		superblock->root_group_object_header_address =
			readValue(reader, superblock->size_of_offsets);

		if (superblock->base_address != 0) {
			mylog("base address is not null\n");
			return MYSOFA_UNSUPPORTED_FORMAT;
		}

		if (fseek(reader->fhd, 0L, SEEK_END))
			return errno;

		if ((long)superblock->end_of_file_address != ftell(reader->fhd)) {
			mylog("file size mismatch\n");
			return MYSOFA_INVALID_FORMAT;
		}

		/* end of superblock */

		/* seek to first object */
		if (fseek(reader->fhd, superblock->root_group_object_header_address,
			SEEK_SET)) {
			mylog("cannot seek to first object at %" PRId64 "\n",
				superblock->root_group_object_header_address);
			return errno;
		}

		return dataobjectRead(reader, &superblock->dataobject, NULL);
	}

	int superblockRead0or1(struct READER *reader, struct SUPERBLOCK *superblock,
		int version) {
		if (fgetc(reader->fhd) !=
			0) /* Version Number of the File’s Free Space Information */
			return MYSOFA_INVALID_FORMAT;

		if (fgetc(reader->fhd) !=
			0) /* Version Number of the Root Group Symbol Table Entry */
			return MYSOFA_INVALID_FORMAT;

		if (fgetc(reader->fhd) != 0)
			return MYSOFA_INVALID_FORMAT;

		if (fgetc(reader->fhd) !=
			0) /* Version Number of the Shared Header Message Format */
			return MYSOFA_INVALID_FORMAT;

		superblock->size_of_offsets = (uint8_t)fgetc(reader->fhd);
		superblock->size_of_lengths = (uint8_t)fgetc(reader->fhd);

		if (fgetc(reader->fhd) != 0)
			return MYSOFA_INVALID_FORMAT;

		if (superblock->size_of_offsets < 2 || superblock->size_of_offsets > 8 ||
			superblock->size_of_lengths < 2 || superblock->size_of_lengths > 8) {
			mylog("size of offsets and length is invalid: %d %d\n",
				superblock->size_of_offsets, superblock->size_of_lengths);
			return MYSOFA_UNSUPPORTED_FORMAT;
		}

		/**
		00000000  89 48 44 46 0d 0a 1a 0a  02 08 08 00 00 00 00 00  |.HDF............|
		00000010  00 00 00 00 ff ff ff ff  ff ff ff ff 72 4d 02 00  |............rM..|
		00000020  00 00 00 00 30 00 00 00  00 00 00 00 e5 a2 a2 a2  |....0...........|
		00000000  89 48 44 46 0d 0a 1a 0a  00 00 00 00 00 08 08 00  |.HDF............|
		00000010  04 00 10 00 01 00 00 00  00 00 00 00 00 00 00 00  |................|
		00000020  ff ff ff ff ff ff ff ff  36 41 9c 00 00 00 00 00  |........6A......|
		00000030  ff ff ff ff ff ff ff ff  00 00 00 00 00 00 00 00  |................|
		00000040  60 00 00 00 00 00 00 00  01 00 00 00 00 00 00 00  |`...............|
		00000050  88 00 00 00 00 00 00 00  a8 02 00 00 00 00 00 00  |................|
		00000060  01 00 1c 00 01 00 00 00  18 00 00 00 00 00 00 00  |................|
		00000070  10 00 10 00 00 00 00 00  20 03 00 00 00 00 00 00  |........ .......|
		00000080  c8 07 00 00 00 00 00 00  54 52 45 45 00 00 02 00  |........TREE....|
		00000000  89 48 44 46 0d 0a 1a 0a  00 00 00 00 00 08 08 00  |.HDF............|
		00000010  04 00 10 00 00 00 00 00  00 00 00 00 00 00 00 00  |................|
		00000020  ff ff ff ff ff ff ff ff  a6 e6 11 00 00 00 00 00  |................|
		00000030  ff ff ff ff ff ff ff ff  00 00 00 00 00 00 00 00  |................|
		00000040  60 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  |`...............|
		00000050  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  |................|
		00000060  4f 48 44 52 02 0d 45 02  02 22 00 00 00 00 00 03  |OHDR..E.."......|
		*/

		/*int groupLeafNodeK = */ readValue(reader, 2);
		/*int groupInternalNodeK = */ readValue(reader, 2);

		if (readValue(reader, 4) != 0) { // has been 1 not 0 in Steinberg file
			mylog("File Consistency Flags are not zero\n");
			return MYSOFA_UNSUPPORTED_FORMAT;
		}

		if (version == 1) {
			readValue(reader, 2); /* Indexed Storage Internal Node K */
			readValue(reader, 2); /* Reserved (zero) */
		}

		superblock->base_address = readValue(reader, superblock->size_of_offsets);
		if (superblock->base_address != 0) {
			mylog("base address is not null\n");
			return MYSOFA_UNSUPPORTED_FORMAT;
		}

		readValue(reader,
			superblock->size_of_offsets); /* Address of File Free space Info */

		superblock->end_of_file_address =
			readValue(reader, superblock->size_of_offsets);

		readValue(reader,
			superblock->size_of_offsets); /* Driver Information Block Address */

		readValue(reader, superblock->size_of_offsets); /* Link Name Offset */

		superblock->root_group_object_header_address = readValue(
			reader, superblock->size_of_offsets); /* Object Header Address */

		uint64_t cache_type = readValue(reader, 4);
		switch (cache_type) {
		case 0:
		case 1:
		case 2:
			break;
		default:
			mylog("cache type must be 0,1, or 2 not %" PRIu64 "\n", cache_type);
			return MYSOFA_UNSUPPORTED_FORMAT;
		}

		if (fseek(reader->fhd, 0L, SEEK_END))
			return errno;

		if ((long)superblock->end_of_file_address != ftell(reader->fhd)) {
			mylog("file size mismatch\n");
		}
		/* end of superblock */

		/* seek to first object */
		if (fseek(reader->fhd, superblock->root_group_object_header_address,
			SEEK_SET)) {
			mylog("cannot seek to first object at %" PRId64 "\n",
				superblock->root_group_object_header_address);
			return errno;
		}

		return dataobjectRead(reader, &superblock->dataobject, NULL);
	}

	int superblockRead(struct READER *reader, struct SUPERBLOCK *superblock) {
		char buf[8];
		memset(superblock, 0, sizeof(*superblock));

		/* signature */
		if (fread(buf, 1, 8, reader->fhd) != 8 ||
			strncmp("\211HDF\r\n\032\n", buf, 8)) {
			mylog("file does not have correct signature");
			return MYSOFA_INVALID_FORMAT;
		}

		/* read version of superblock, must be 0,1,2, or 3 */
		int version = fgetc(reader->fhd);

		switch (version) {
		case 0:
		case 1:
			return superblockRead0or1(reader, superblock, version);
		case 2:
		case 3:
			return superblockRead2or3(reader, superblock);
		default:
			mylog("superblock must have version 0, 1, 2, or 3 but has %d\n", version);
			return MYSOFA_INVALID_FORMAT;
		}
	}

	void superblockFree(struct READER *reader, struct SUPERBLOCK *superblock) {
		dataobjectFree(reader, &superblock->dataobject);
	}

#else
extern int to_avoid_iso_compiler_warning_when_there_are_no_symbols;
#endif /* SAF_ENABLE_SOFA_READER_MODULE */
